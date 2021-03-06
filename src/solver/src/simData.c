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
#include "solver/simData.h"
#include "solver/quadData.h"
#include "solver/refine.h"
#include "solver/coarsen.h"
#include "solver/gradients.h"
#include "solver/paramfile.h"
#include "aux/dbg.h"

#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#endif


/***********************************************************
* init_simData()
*-----------------------------------------------------------
* Initializes the simulation data structure
***********************************************************/
SimData_t *init_simData(int          argc, 
                        char        *argv[], 
                        octInitFun   usrInitFun,
                        octRefineFun usrRefineFun,
                        octCoarseFun usrCoarseFun)
{
  /*--------------------------------------------------------
  | Check that parameter file exists
  --------------------------------------------------------*/
  check(argc > 1,
      "\n\nNo parameter file provided!\n\n  OctFS <Parameter file>\n");
  const char *paramFilePath = argv[1];

  /*--------------------------------------------------------
  | Create data structures
  --------------------------------------------------------*/
  SimData_t *simData = malloc(sizeof(SimData_t));

  simData->simParam    = NULL;
  simData->solverParam = NULL;
  simData->mpiParam    = NULL;
  simData->conn        = NULL;
  simData->p4est       = NULL;

  /*--------------------------------------------------------
  | Init parameter structures 
  --------------------------------------------------------*/
  simData->simParam    = init_simParam(usrInitFun, 
                                       usrRefineFun,
                                       usrCoarseFun);
  simData->solverParam = init_solverParam();
  simData->mpiParam    = init_mpiParam(argc, argv);

  /*--------------------------------------------------------
  | Read parameter file
  --------------------------------------------------------*/
  octParam_readParamfile(simData, paramFilePath);

  /*--------------------------------------------------------
  | Load p4est mesh connectivity
  --------------------------------------------------------*/
#ifndef P4_TO_P8
  simData->conn = p4est_connectivity_new_periodic();
#else
  simData->conn = p8est_connectivity_new_periodic();
#endif

  /*--------------------------------------------------------
  | create p4est structure
  --------------------------------------------------------*/
  SolverParam_t *solverParam = simData->solverParam;

  p4est_init(NULL, SC_LP_PRODUCTION);
  P4EST_GLOBAL_PRODUCTIONF(
      "\n\nOctFS - Octree based flow solver. Compiled for %dD.\n\n",
      P4EST_DIM);
  simData->p4est = p4est_new_ext(simData->mpiParam->mpiComm,
                                 simData->conn,
                                 solverParam->nQuadMPU,
                                 solverParam->minRefLvl,
                                 solverParam->fillUniform,
                                 sizeof(QuadData_t),
                                 init_quadData,
                                 (void *) (simData));

  /*--------------------------------------------------------
  | Estimate global mesh attributes
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, NULL, NULL,
                estimateMeshAttributes, // cell callback
                NULL,               // face callback
#ifdef P4_TO_P8
                NULL,               // edge callback
#endif
                NULL);              // corner callback*/

  /*--------------------------------------------------------
  | Exchange data global among processes
  --------------------------------------------------------*/
  exchangeGlobMeshAttrib(simData);

  /*--------------------------------------------------------
  | Init p4est ghost data structure
  --------------------------------------------------------*/
  simData->ghost = p4est_ghost_new(simData->p4est, 
                                   P4EST_CONNECT_FULL);
  simData->ghostData = P4EST_ALLOC(QuadData_t, 
                          simData->ghost->ghosts.elem_count);
  p4est_ghost_exchange_data(simData->p4est, 
                            simData->ghost, 
                            simData->ghostData);

  /*--------------------------------------------------------
  | Initial calculation of gradients
  --------------------------------------------------------*/
  int idx;
  for (idx = 0; idx < OCT_MAX_VARS; idx++)
    computeGradients(simData, idx); 

  if (solverParam->adaptGrid == TRUE)
  {
    /*------------------------------------------------------
    | Initial refinement 
    ------------------------------------------------------*/
    p4est_refine(simData->p4est,
                 solverParam->recursive,
                 globalRefinement,
                 init_quadData);

    p4est_coarsen(simData->p4est,
                  solverParam->recursive,
                  globalCoarsening,
                  init_quadData);

    p4est_balance(simData->p4est,
                  P4EST_CONNECT_FACE,
                  init_quadData);

    p4est_ghost_destroy(simData->ghost);
    P4EST_FREE(simData->ghostData);
    simData->ghost = NULL;
    simData->ghostData = NULL;
  }

  simData->ghost = p4est_ghost_new(simData->p4est, 
                                   P4EST_CONNECT_FULL);
  simData->ghostData = P4EST_ALLOC(QuadData_t, 
                          simData->ghost->ghosts.elem_count);
  p4est_ghost_exchange_data(simData->p4est, 
                            simData->ghost, 
                            simData->ghostData);

  for (idx = 0; idx < OCT_MAX_VARS; idx++)
    computeGradients(simData, idx); 

  return simData;

error:
  return NULL;

} /* init_simData() */

/***********************************************************
* init_simParam()
*-----------------------------------------------------------
* Initializes the simulation parameter structure
***********************************************************/
SimParam_t *init_simParam(octInitFun   usrInitFun,
                          octRefineFun usrRefineFun,
                          octCoarseFun usrCoarseFun)
{
  int i;

  SimParam_t *simParam = malloc(sizeof(SimParam_t));

  simParam->volume_glob     = 0.0;
  simParam->volume_loc      = 0.0;

  simParam->timestep      = 5e-3;
  simParam->simTimeTot    = 1.0; //1.0; //5e-3;
  simParam->simTime       = 0.0;

  simParam->tempScheme     = CRANK_NICOLSON;
  simParam->tempFluxFac[0] = 0.0;
  simParam->tempFluxFac[1] = 1.0;
  simParam->tempFluxFac[2] = 0.5;
  

  simParam->viscosity     = 1e-5;
  simParam->ref_length    = 1.0;
  simParam->ref_velocity  = 1.0;
  simParam->ref_pressure  = 0.0;

  simParam->usrInitFun   = usrInitFun;
  simParam->usrRefineFun = usrRefineFun;
  simParam->usrCoarseFun = usrCoarseFun;

  /*--------------------------------------------------------
  | Temporary values
  --------------------------------------------------------*/
  simParam->tmp_xId  = -1;
  simParam->tmp_AxId = -1;

  simParam->tmp_fluxFac = 0.0;

  simParam->tmp_sbufVec0 = -1;
  simParam->tmp_sbufVec1 = -1;
  simParam->tmp_sbufProd = -1;

  for (i = 0; i < PARAM_BUF_VARS; i++)
    simParam->sbuf[i] = 0.0;

  return simParam;

} /* init_simParam() */

/***********************************************************
* init_solverParam()
*-----------------------------------------------------------
* Initializes the solver parameter structure
***********************************************************/
SolverParam_t *init_solverParam(void)
{
  SolverParam_t *solverParam = malloc(sizeof(SolverParam_t));

  solverParam->epsilon = 1.0E-06;


  // Path to export directory
  solverParam->io_exportDir = "./";
  // Prefix for export files
  solverParam->io_exportPrefix = "TestRun";


  // Number of quadrants per MPU
  solverParam->nQuadMPU = 0;
  // Minimum level of refinement for initialization
  solverParam->minRefLvl = 5;
  // Maximum level of refinement
  solverParam->maxRefLvl = 6;
  // Fill uniform
  solverParam->fillUniform = TRUE;
  // Recursive refinement
  solverParam->recursive = TRUE;
  // Re-Partition on coarsening
  solverParam->partForCoarsen = TRUE;
  // Turn on/off automatic grid adaptation
  solverParam->adaptGrid = TRUE;

  // Global refinement error for passive scalar 
  solverParam->refErr_scalar    = 0.05;
  // Global refinement error for pressure  
  solverParam->refErr_pressure  = 1.0E-03;


  // Number of timesteps between refinement periods
  solverParam->refinePeriod = 1;
  // Numer of timesteps between repartitioning
  solverParam->repartitionPeriod = 1;
  // Number of timesteps between solution writes
  solverParam->writePeriod = 10;

  return solverParam;

} /* init_solverParam() */

/***********************************************************
* init_mpiParam()
*-----------------------------------------------------------
* Initializes the solver parameter structure
***********************************************************/
MPIParam_t *init_mpiParam(int argc, char *argv[])
{
  MPIParam_t *mpiParam = NULL;
  mpiParam             = malloc(sizeof(MPIParam_t));

  int mpi_return;
  
  mpi_return = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpi_return);

  mpiParam->mpiComm = sc_MPI_COMM_WORLD;

  return mpiParam;

} /* init_mpiParam() */




/***********************************************************
* destroy_simData()
*-----------------------------------------------------------
* Frees all memory of a SimData structure
***********************************************************/
void destroy_simData(SimData_t *simData)
{

  destroy_simParam(simData->simParam);
  destroy_solverParam(simData->solverParam);

  if (simData->p4est != NULL)
    p4est_destroy(simData->p4est);

  if (simData->conn != NULL)
    p4est_connectivity_destroy(simData->conn);

  if (simData->mpiParam != NULL)
  destroy_mpiParam(simData->mpiParam);

  free(simData);


} /* destroy_simData() */

/***********************************************************
* destroy_simParam()
*-----------------------------------------------------------
* Frees all memory of a SimParam structure
***********************************************************/
void destroy_simParam(SimParam_t *simParam)
{
  free(simParam);

} /* destroy_simParam() */

/***********************************************************
* destroy_solverParam()
*-----------------------------------------------------------
* Frees all memory of a SolverParam structure
***********************************************************/
void destroy_solverParam(SolverParam_t *solverParam)
{
  free (solverParam);

} /* destroy_solverParam() */

/***********************************************************
* destroy_mpiParam()
*-----------------------------------------------------------
* Frees all memory of a MPIParam structure
***********************************************************/
void destroy_mpiParam(MPIParam_t *mpiParam)
{
  int mpi_return = sc_MPI_Finalize();
  SC_CHECK_MPI(mpi_return);

  free(mpiParam);
  
} /* destroy_mpiParam() */

/***********************************************************
* estimateMeshAttributes()
*-----------------------------------------------------------
* Estimate global mesh attributes such as number of elements
* e.g.
*
*   -> p4est_iter_volume_t callback function
***********************************************************/
void estimateMeshAttributes(p4est_iter_volume_info_t *info,
                            void *user_data)
{
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;
  SimData_t  *simData  = (SimData_t*)info->p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  simParam->volume_loc += quadData->volume;

} /* estimateMeshAttributes() */

/***********************************************************
* exchangeGlobMeshAttrib()
*-----------------------------------------------------------
* Function to exchange global mesh attributes among all
* processes
***********************************************************/
void exchangeGlobMeshAttrib(SimData_t *simData)
{
  sc_MPI_Allreduce(&simData->simParam->volume_loc,
                   &simData->simParam->volume_glob,
                   1,
                   sc_MPI_DOUBLE,
                   sc_MPI_SUM,
                   simData->mpiParam->mpiComm);

} /* exchangeGlobMeshAttrib() */
