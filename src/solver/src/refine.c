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
* Static variables
***********************************************************/
static int varIdx;

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
octDouble calcSqrErr_scalar(p4est_quadrant_t *q)
{
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;
  QuadFlowData_t *flowData = &quadData->flowData;
  QuadGeomData_t *geomData = &quadData->geomData;
  octDouble *grad_s = flowData->grad_vars[varIdx];
  octDouble  vol    = geomData->volume;

  int i;
  octDouble diff2 = 0.;
  /* use the approximate derivative to estimate the L2 error */
  for (i = 0; i < P4EST_DIM; i++) {
    diff2 += grad_s[i] * grad_s[i] * (1. / 12.) * vol;
  }

  return diff2;
}


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
  QuadGeomData_t *geomData = &quadData->geomData;

  SimData_t *simData = (SimData_t *) p4est->user_pointer;

  octDouble globErr   = simData->solverParam->refErr_scalar;
  octDouble globErr2  = globErr * globErr;

  octDouble vol = geomData->volume; 

  varIdx = IS;
  octDouble err2 = calcSqrErr_scalar(q);

  if (err2 > globErr2 * vol) {
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
