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
* init_quadData()
*-----------------------------------------------------------
* Initializes the quad data structure
***********************************************************/
void init_quadData(p4est_t *p4est,
                   p4est_topidx_t which_tree,
                   p4est_quadrant_t *q)
{
  SimData_t  *simData  = (SimData_t *) p4est->user_pointer;
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;

  QuadGeomData_t *geomData = &(quadData->geomData);
  QuadFlowData_t *flowData = &(quadData->flowData);

#ifdef P4_TO_P8
  init_quadGeomData3d(p4est, which_tree, q, geomData);
#else
  init_quadGeomData2d(p4est, which_tree, q, geomData);
#endif


  init_quadFlowData(p4est, which_tree, q, flowData);

  /*--------------------------------------------------------
  | Apply user-defined initialization function 
  --------------------------------------------------------*/
  if (simData->simParam->usrInitFun != NULL)
  {
    simData->simParam->usrInitFun(quadData);
  }


} /* init_quadData() */


/***********************************************************
* init_quadFlowData()
*-----------------------------------------------------------
* Initializes the quad flow data structure
***********************************************************/
void init_quadFlowData(p4est_t *p4est,
                       p4est_topidx_t which_tree,
                       p4est_quadrant_t *q, 
                       QuadFlowData_t *flowData)
{
  int i,j;

  for (i = 0; i < OCT_MAX_VARS; i++)
  {
    flowData->vars[i]     = 0.0;
    flowData->vars_buf[i] = 0.0;

    for (j = 0; j < P4EST_DIM; j++)
    {
      flowData->grad_vars[i][j] = 0.0;
    }
  }

} /* init_quadFlowData() */


/***********************************************************
* init_quadGeomData()
*-----------------------------------------------------------
* Initializes the quad geometry data structure
*
*                 n[3]
*        V[2]<-------------V[3]
*         |                 ^
*         |                 |
*         |                 |
*     n[2]|                 | n[1]
*         |                 |
*   y     |                 |
*   |     v                 |
*   |    V[0]------------->V[1]
*   |              n[0]
*   ------>x
*  
*
***********************************************************/
void init_quadGeomData2d(p4est_t *p4est,
                         p4est_topidx_t which_tree,
                         p4est_quadrant_t *q, 
                         QuadGeomData_t *geomData)
{
  int i,j;

  /*--------------------------------------------------------
  | Get 2D vertex coordinates 
  --------------------------------------------------------*/
  p4est_qcoord_t length = P4EST_QUADRANT_LEN(q->level);

  octDouble (*xyz)[2] = geomData->xyz;

  // V[0]
  p4est_qcoord_to_vertex(p4est->connectivity,
                         which_tree,
                         q->x,
                         q->y,
                         xyz[0]);

  // V[1]
  p4est_qcoord_to_vertex(p4est->connectivity,
                         which_tree,
                         q->x + length,
                         q->y,
                         xyz[1]);

  // V[2]
  p4est_qcoord_to_vertex(p4est->connectivity,
                         which_tree,
                         q->x,
                         q->y + length,
                         xyz[2]);

  // V[3]
  p4est_qcoord_to_vertex(p4est->connectivity,
                         which_tree,
                         q->x + length,
                         q->y + length,
                         xyz[3]);

  /*--------------------------------------------------------
  | Compute quadrant area
  --------------------------------------------------------*/
  geomData->volume = ( xyz[0][0] * xyz[1][1]
                     + xyz[1][0] * xyz[2][1]
                     + xyz[2][0] * xyz[3][1]
                     + xyz[3][0] * xyz[0][1] )
                    -( xyz[1][0] * xyz[0][1]
                     + xyz[2][0] * xyz[1][1]
                     + xyz[3][0] * xyz[2][1]
                     + xyz[0][0] * xyz[3][1] ) * 0.5;

  /*--------------------------------------------------------
  | Compute quadrant centroid
  --------------------------------------------------------*/
  geomData->centroid[0] = 0.25 * ( xyz[0][0] + xyz[1][0]
                                 + xyz[2][0] + xyz[3][0] );
  geomData->centroid[1] = 0.25 * ( xyz[0][1] + xyz[1][1]
                                 + xyz[2][1] + xyz[3][1] );

  /*--------------------------------------------------------
  | Compute quadrant normals
  --------------------------------------------------------*/
  geomData->normals[0][0] = -xyz[1][1] + xyz[0][1];
  geomData->normals[0][1] =  xyz[1][0] - xyz[0][0];

  geomData->normals[1][0] = -xyz[3][1] + xyz[1][1];
  geomData->normals[1][1] =  xyz[3][0] - xyz[1][0];

  geomData->normals[2][0] = -xyz[0][1] + xyz[2][1];
  geomData->normals[2][1] =  xyz[0][0] - xyz[2][0];

  geomData->normals[3][0] = -xyz[2][1] + xyz[3][1];
  geomData->normals[3][1] =  xyz[2][0] - xyz[3][0];

  /*--------------------------------------------------------
  | Compute quadrant normal centroids
  --------------------------------------------------------*/
  geomData->face_centroids[0][0] = 0.5*(xyz[1][0]+xyz[0][0]);
  geomData->face_centroids[0][1] = 0.5*(xyz[1][1]+xyz[0][1]);

  geomData->face_centroids[1][0] = 0.5*(xyz[3][0]+xyz[1][0]);
  geomData->face_centroids[1][1] = 0.5*(xyz[3][1]+xyz[1][1]);

  geomData->face_centroids[2][0] = 0.5*(xyz[0][0]+xyz[2][0]);
  geomData->face_centroids[2][1] = 0.5*(xyz[0][1]+xyz[2][1]);

  geomData->face_centroids[3][0] = 0.5*(xyz[2][0]+xyz[3][0]);
  geomData->face_centroids[3][1] = 0.5*(xyz[2][1]+xyz[3][1]);


} /* init_quadGeomData() */

