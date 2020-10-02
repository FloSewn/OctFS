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
#include "solver/util.h"
#include "solver/coarsen.h"
#include "solver/quadData.h"
#include "solver/simData.h"

/***********************************************************
* calcSqrErr()
*-----------------------------------------------------------
* This function estimates the approximation error 
* on a quadrant for a specified variable.
*
* The estimate is computed by integrating the difference of a 
* constant approximation at the midpoint and a linear 
* approximation that interpolates at the midpoint.
*
* The square of the error estiamte for the state variables 
* contained is returned.
*
***********************************************************/
octDouble calcSqrErr(p4est_quadrant_t *q, int varIdx)
{
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;

  octDouble *grad_s = quadData->grad_vars[varIdx];

  octDouble l = (octDouble) P4EST_ROOT_LEN;
  octDouble h = (octDouble) P4EST_QUADRANT_LEN(q->level) / l;
  
  /*--------------------------------------------------------
  | use approximate derivative to estimate L2 error of 
  | associated energy
  --------------------------------------------------------*/
  int i;
  octDouble k = 0.0;
  octDouble v = 1.0;
  
  for (i = 0; i < P4EST_DIM; i++) 
  {
    k += 0.5 * grad_s[i] * h;
    v *= h;
  }

  return v * (k * k * k * k);


} /* calcSqrErr() */


/***********************************************************
* refinement_scalarError()
*-----------------------------------------------------------
* Function to calculate the error estimate for the mesh
* refinement 
***********************************************************/
int refinement_scalarError(p4est_t          *p4est, 
                           p4est_topidx_t    which_tree,
                           p4est_quadrant_t *q)
{
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;

  SimData_t *simData  = (SimData_t *) p4est->user_pointer;
  octDouble  globErr  = simData->solverParam->refErr_scalar;
  octDouble  globErr2 = globErr * globErr;

  octDouble vol = quadData->volume; 

  octDouble err2 = calcSqrErr(q, IS);

  if (err2 > globErr2 * globErr2 * vol) 
  {
    return 1;
  }

  return 0;

} /* refinement_scalarError() */


/***********************************************************
* globalRefinement()
*-----------------------------------------------------------
* This is the general function to control the refinement
* of trees.
***********************************************************/
int globalRefinement(p4est_t          *p4est,
                     p4est_topidx_t    which_tree,
                     p4est_quadrant_t *q)
{
  SimData_t  *simData  = (SimData_t*) p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  int refine = 0;

  refine |= refinement_scalarError(p4est, which_tree, q);

  if (simParam->usrRefineFun != NULL)
  {
    refine |= simParam->usrRefineFun(p4est, which_tree, q);
  }

  return refine;

} /* globalRefinement() */
