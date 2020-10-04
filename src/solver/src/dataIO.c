/*
* This file is part of OctFS. 
* OctFS is a finite-volume flow solver with adaptive
* mesh refinement written in C, which is based on 
* the p4est library.
*
* Copyright (C) 2020 Florian Setzwein 
*
* OctFS is free software; you can redistribute it and/or 
* modify it under the terms of the GNU General Public 
* License as published by the Free Software Foundation; 
* either version 2 of the License, or (at your option) 
* any later version.
*
* OctFS is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied 
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
* PURPOSE.  See the GNU General Public License for more 
* details.
*
* You should have received a copy of the GNU General 
* Public License along with OctFS; if not, write to the 
* Free Software Foundation, Inc., 51 Franklin Street, 
* Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include "solver/typedefs.h"
#include "solver/simData.h"
#include "solver/quadData.h"
#include "solver/dataIO.h"
#include "solver/util.h"
#include "aux/dbg.h"

/***********************************************************
* Solver variable name definition
***********************************************************/
static char varNames[OCT_MAX_VARS][OCT_VARNAME_LENGTH] = 
{
  "solver_Ax",
  "solver_b",
  "solver_vn",
  "solver_r",
  "solver_r0",
  "solver_p",
  "solver_v",
  "solver_h",
  "solver_s",
  "solver_t",
  "solver_res",
  "density",
  "x_velocity",
  "y_velocity",
  "z_velocity",
  "pressure",
  "passive_scalar"
}; 

/***********************************************************
* Variable index
***********************************************************/
static int io_idx = -1;

/***********************************************************
* Function to concatenate two strings.
* Appends <str_2> to the end of <str_1>.
*-----------------------------------------------------------
* Arguments:
* *str_1, str_2 : strings to concatenate
*
***********************************************************/
char *concat_string(const char *str_1, const char *str_2)
{
  const size_t str_1_len  = strlen(str_1);
  const size_t str_2_len  = strlen(str_2);
  const size_t total_len = str_1_len + str_2_len;

  char *const strBuf = malloc(total_len + 1);
  check_mem(strBuf);

  strcpy(strBuf, str_1);
  strcpy(strBuf + str_1_len, str_2);

  return strBuf;

error:
  return NULL;

}

/***********************************************************
* Callback function for interpolating the solution from 
* quadrant midpoints to corners.
*
* The function p4est_iterate() takes as an argument a 
* p4est_iter_volume_t callback function, which it executes at 
* every local quadrant (see p4est_iterate.h).  
* This function matches the p4est_iter_volume_t prototype.
*
* Use the callback function to interpolate the state
* variable to the corners, and write those corners into an 
* array so that they can be written out.
*-----------------------------------------------------------
* Arguments:
* *info       : the information about this quadrant that 
*               has been  populated by p4est_iterate()
* user_data   : the user_data that was given as an argument to
*               p4est_iterate: in this case, it points to the
*               array of corner values that we want to write.
*               The values for the corner of the quadrant
*               described by a info are written during the
*               execution of the callback.
*
***********************************************************/
static void interpSolution(p4est_iter_volume_info_t *info, 
                           void *user_data)
{
  /*--------------------------------------------------------
  | we passed the array of values to fill as the user_data 
  | in the call to p4est_iterate 
  |-------------------------------------------------------*/
  sc_array_t         *var_interp = (sc_array_t *) user_data;      
  p4est_t            *p4est      = info->p4est;
  p4est_quadrant_t   *q          = info->quad;
  p4est_topidx_t      which_tree = info->treeid;

  /*--------------------------------------------------------
  | local_id is the index of q *within its tree's numbering.  
  | We want to convert it its index for all the quadrants on 
  | this process, which we do below.
  |-------------------------------------------------------*/
  p4est_locidx_t  local_id = info->quadid;  
  QuadData_t     *quadData = (QuadData_t *) q->p.user_data;
  p4est_tree_t   *tree     = p4est_tree_array_index(p4est->trees, 
                                                    which_tree);

  // now the id is relative to the MPI process 
  local_id += tree->quadrants_offset;   
  
  // each local quadrant has 2^d (P4EST_CHILDREN) 
  // values in var_interp 
  p4est_locidx_t arrayoffset = P4EST_CHILDREN * local_id;      
  //double h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;

  double  this_u;
  double *this_u_ptr;
  int     i, j;

  for (i = 0; i < P4EST_CHILDREN; i++) 
  {
    this_u = quadData->vars[io_idx];
    //this_u = quadData->b[io_idx] / quadData->volume;
    //this_u = quadData->grad_vars[io_idx][1];

    /*------------------------------------------------------
    | loop over the derivative components and 
    | linearly interpolate from the midpoint to the corners 
    |-----------------------------------------------------*/
    for (j = 0; j < P4EST_DIM; j++) 
    {
      /*----------------------------------------------------
      | In order to know whether the direction from the 
      | midpoint to the corner is negative or positive, we 
      | take advantage of the fact that the corners are 
      | in z-order.  
      | If i is an odd number, it is 
      | on the +x side; if it is even, it is on the -x side.  
      | If (i / 2) is an odd number, it is on the +y sidea.. 
      |---------------------------------------------------*/
      // this_u += (h / 2) * data->grad_u[j] * ((i & (1 << j)) ? 1. : -1.);
    }

    this_u_ptr = (double *) sc_array_index(var_interp, 
                                           arrayoffset + i);
    this_u_ptr[0] = this_u;
  }

} /* interpSolution(...) */


/***********************************************************
* writeSolutionVtk()
*-----------------------------------------------------------
* Function to write the solution of a single timestep 
* to vtk format
***********************************************************/
void writeSolutionVtk(SimData_t *simData, int step)
{
  p4est_t       *p4est       = simData->p4est;
  SolverParam_t *solverParam = simData->solverParam;

  char                 filename[BUFSIZ] = "";
  int                  retval;
  p4est_vtk_context_t *context;

  snprintf(filename, BUFSIZ, "%s_%04d", 
      concat_string(solverParam->io_exportDir,
                    solverParam->io_exportPrefix), step);
  P4EST_GLOBAL_PRODUCTIONF("Writing results file: %s\n", 
      filename);

  /*--------------------------------------------------------
  | create a vector with one value for the corner of 
  | every local quadrant (the number of children is always 
  | the same as the number of corners)
  |-------------------------------------------------------*/
  p4est_locidx_t numquads = p4est->local_num_quadrants;
  int nEntries            = numquads * P4EST_CHILDREN;
  sc_array_t *var_interp[OCT_MAX_VARS];
  
  for (io_idx = 0; io_idx < OCT_MAX_VARS; io_idx++)
    var_interp[io_idx] = sc_array_new_size(sizeof (double), 
                                           nEntries);


  /*--------------------------------------------------------
  | create VTK output context and set its parameters
  |-------------------------------------------------------*/
  context = p4est_vtk_context_new(p4est, filename);

  /* quadrant at scale 1.0 */
  p4est_vtk_context_set_scale(context, 0.95);  

  /* begin writing the output files */
  context = p4est_vtk_write_header(context);

  SC_CHECK_ABORT(context != NULL,
                 P4EST_STRING "_vtk: Error writing vtk header");

  /* do not write the tree id's of each quadrant */
  context = p4est_vtk_write_cell_dataf(context, 0, 
                                       1, // write refinelevel 
                                       1, // write mpi id      
                                       0, // wrap mpi rank     
                                       0, // custom scalar data
                                       0, // custom vector data
                                       context);      

  SC_CHECK_ABORT(context != NULL,
                 P4EST_STRING "_vtk: Error writing cell data");

  /*--------------------------------------------------------
  | Interpolate field data to points
  |-------------------------------------------------------*/
  for (io_idx = 0; io_idx < OCT_MAX_VARS; io_idx++)
  {
    p4est_iterate(p4est, 
                  NULL,   
                  (void *) var_interp[io_idx],     
                  interpSolution,
                  NULL,          
#ifdef P4_TO_P8
                  NULL,         
#endif
                  NULL);         
  }

  /*--------------------------------------------------------
  | Set field data for output file
  |-------------------------------------------------------*/
  context = p4est_vtk_write_point_dataf(context, 
#ifdef P4_TO_P8
                              OCT_MAX_VARS - OCT_SOLVER_VARS,
#else
                              OCT_MAX_VARS - OCT_SOLVER_VARS - 1,
#endif
                              0, 
                              varNames[0], var_interp[0], 
                              varNames[2], var_interp[2], 
#ifdef P4_TO_P8
                              varNames[3], var_interp[3], 
#endif
                              varNames[3], var_interp[3], 
                              varNames[5], var_interp[5], 
                              varNames[7], var_interp[7], 
                              context);    
  SC_CHECK_ABORT(context != NULL,
                 P4EST_STRING "_vtk: Error writing field data");

  /*--------------------------------------------------------
  | Finalize vtk file
  |-------------------------------------------------------*/
  retval = p4est_vtk_write_footer(context);

  SC_CHECK_ABORT(!retval, P4EST_STRING "_vtk: Error writing footer");

  for (io_idx = 0; io_idx < OCT_MAX_VARS; io_idx++)
    sc_array_destroy(var_interp[io_idx]);


} /* writeSolutionVtk() */
