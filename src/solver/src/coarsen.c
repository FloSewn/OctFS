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
#include "solver/refine.h"
#include "solver/quadData.h"
#include "solver/simData.h"

/***********************************************************
* coarsening_scalarError()
*-----------------------------------------------------------
* Function to calculate the error estimate for the mesh
* coarsening 
***********************************************************/
int coarsening_scalarError(p4est_t *p4est,
                           p4est_topidx_t which_tree,
                           p4est_quadrant_t * children[])
{
  SimData_t *simData  = (SimData_t *) p4est->user_pointer;
  octDouble  globErr  = simData->solverParam->refErr_scalar;
  octDouble  globErr2 = globErr * globErr;

  /*--------------------------------------------------------
  | Compute average of children values
  ---------------------------------------------------------*/
  QuadData_t     *quadData;
  QuadGeomData_t *geomData;
  QuadFlowData_t *flowData;

  int i;
  octDouble var_p = 0.;
  octDouble vol_p = 0.;

  for (i = 0; i < P4EST_CHILDREN; i++) 
  {
    quadData = (QuadData_t *) children[i]->p.user_data;
    geomData = &quadData->geomData;
    flowData = &quadData->flowData;

    const octDouble vol = geomData->volume;
    const octDouble var = flowData->vars[IS];

    var_p += vol * var;
    vol_p += vol;
  }

  /*--------------------------------------------------------
  | Compute error
  ---------------------------------------------------------*/
  octDouble err2 = 0.0;

  for (i = 0; i < P4EST_CHILDREN; i++) 
  {
    quadData = (QuadData_t *) children[i]->p.user_data;
    geomData = &quadData->geomData;
    flowData = &quadData->flowData;

    const octDouble err2_c = calcSqrErr(children[i], IS);

    const octDouble vol = geomData->volume;

    if (err2_c > globErr2 * vol)
      return 0;

    err2 += err2_c;

    const octDouble var = flowData->vars[IS];

    const octDouble diff = var_p - var;
    const octDouble diff2 = diff * diff;

    err2 += vol * diff2;
  }

  if (err2 < globErr2 * vol_p)
    return 1;

  return 0;

} /* coarsening_scalarError() */

/***********************************************************
* globalCoarsening()
*-----------------------------------------------------------
* This is the general function to control the coarsening
* of trees.
***********************************************************/
int globalCoarsening(p4est_t          *p4est,
                     p4est_topidx_t    which_tree,
                     p4est_quadrant_t *children[])
{
  SimData_t  *simData  = (SimData_t*) p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  int coarsen = 0;

  coarsen |= coarsening_scalarError(p4est, 
                                    which_tree, 
                                    children);


  if (simParam->usrCoarseFun != NULL)
  {
    coarsen |= simParam->usrCoarseFun(p4est, 
                                      which_tree, 
                                      children);
  }

  return coarsen;

} /* globalCoarsening() */
