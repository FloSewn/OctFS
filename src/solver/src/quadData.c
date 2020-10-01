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

#ifdef P4_TO_P8
  init_quadGeomData3d(p4est, which_tree, q, quadData);
#else
  init_quadGeomData2d(p4est, which_tree, q, quadData);
#endif

  init_quadFlowData(quadData);

  /*--------------------------------------------------------
  | Apply user-defined initialization function as 
  | initialization. Otherwise, interpolate solution
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
void init_quadFlowData(QuadData_t *quadData)
{
  int i,j;

  for (i = 0; i < OCT_MAX_VARS; i++)
  {
    quadData->vars[i]     = 0.0;
    quadData->vars_buf[i] = 0.0;

    for (j = 0; j < P4EST_DIM; j++)
    {
      quadData->grad_vars[i][j] = 0.0;
    }
  }

  /*--------------------------------------------------------
  | buffers for flow solver 
  --------------------------------------------------------*/
  quadData->Ax_p = NULL;

  for (i = 0; i < OCT_MAX_VARS; i++)
  {
    quadData->Ax[i]  = 0.0;
    quadData->b[i]   = 0.0;
    quadData->res[i] = 0.0;
  }


  /*--------------------------------------------------------
  | Set Density to 1
  --------------------------------------------------------*/
  quadData->vars[IRHO] = 1.0;


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
*     n[0]|                 | n[1]
*         |                 |
*   y     |                 |
*   |     v                 |
*   |    V[0]------------->V[1]
*   |              n[2]
*   ------>x
*  
*
***********************************************************/
void init_quadGeomData2d(p4est_t          *p4est,
                         p4est_topidx_t    which_tree,
                         p4est_quadrant_t *q, 
                         QuadData_t       *quadData)
{
  /*--------------------------------------------------------
  | Get 2D vertex coordinates 
  --------------------------------------------------------*/
  p4est_qcoord_t length = P4EST_QUADRANT_LEN(q->level);

  octDouble (*xyz)[2] = quadData->xyz;

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
  octDouble vol = ( xyz[0][0] * xyz[1][1]
                  + xyz[1][0] * xyz[3][1]
                  + xyz[3][0] * xyz[2][1]
                  + xyz[2][0] * xyz[0][1] )
                 -( xyz[0][1] * xyz[1][0]
                  + xyz[1][1] * xyz[3][0]
                  + xyz[3][1] * xyz[2][0]
                  + xyz[2][1] * xyz[0][0] );
  quadData->volume = 0.5 * vol;

  /*--------------------------------------------------------
  | Compute quadrant centroid
  --------------------------------------------------------*/
  quadData->centroid[0] = 0.25 * ( xyz[0][0] + xyz[1][0]
                                 + xyz[2][0] + xyz[3][0] );
  quadData->centroid[1] = 0.25 * ( xyz[0][1] + xyz[1][1]
                                 + xyz[2][1] + xyz[3][1] );

  /*--------------------------------------------------------
  | Compute quadrant normals
  --------------------------------------------------------*/
  quadData->normals[0][0] = xyz[0][1] - xyz[2][1];
  quadData->normals[0][1] = xyz[2][0] - xyz[0][0];

  quadData->normals[1][0] = xyz[3][1] - xyz[1][1];
  quadData->normals[1][1] = xyz[1][0] - xyz[3][0];

  quadData->normals[2][0] = xyz[1][1] - xyz[0][1];
  quadData->normals[2][1] = xyz[0][0] - xyz[1][0];

  quadData->normals[3][0] = xyz[2][1] - xyz[3][1];
  quadData->normals[3][1] = xyz[3][0] - xyz[2][0];

  /*--------------------------------------------------------
  | Compute quadrant normal centroids
  --------------------------------------------------------*/
  quadData->face_centroids[0][0] = 0.5*(xyz[0][0]+xyz[2][0]);
  quadData->face_centroids[0][1] = 0.5*(xyz[0][1]+xyz[2][1]);

  quadData->face_centroids[1][0] = 0.5*(xyz[3][0]+xyz[1][0]);
  quadData->face_centroids[1][1] = 0.5*(xyz[3][1]+xyz[1][1]);

  quadData->face_centroids[2][0] = 0.5*(xyz[1][0]+xyz[0][0]);
  quadData->face_centroids[2][1] = 0.5*(xyz[1][1]+xyz[0][1]);

  quadData->face_centroids[3][0] = 0.5*(xyz[2][0]+xyz[3][0]);
  quadData->face_centroids[3][1] = 0.5*(xyz[2][1]+xyz[3][1]);


} /* init_quadGeomData() */

/***********************************************************
* interpQuadData()
*-----------------------------------------------------------
* Function for initialzing the state variables of incoming
* quadrants from outgoing quadrants. 
*
* The functions p4est_refine_ext(), p4est_coarsen_ext() and
* p4est_balance_ext() take as an arguments a p4est_replace_t
* callback function.
* This function allows to setup the quadrant data of
* incoming quadrants from the data of outgoing quadrants, 
* before the outgoing data is destroyed.
* This function matches the p4est_replace_t prototype.
*
* In this example, we linearly interpolate the state 
* variable of a quadrant that is refined to its children
* and we average the children midpoints that are being 
* coarsened to the paren
*-----------------------------------------------------------
* Arguments:
* *p4est        : the forest
*  which_tree   : the tree in the forest containing a 
*                 child
*  num_outgoing : the number of quadrants that are being
*                 replaced:
*                 Either 1 if a quadrant is being refined 
*                 or P4EST_CHILDREN if a family of 
*                 children are being coarsened
*  outgoing     : the outgoing quadrants
*  num_incoming : the number of quadrants that are being
*                 added:
*                 Either P4EST_CHILDREN if a quadrant is 
*                 being refined, or 1 if a family of 
*                 children are being coarsened.
*  incoming     : quadrants whose data is initialized.
*
***********************************************************/
void interpQuadData(p4est_t          *p4est, 
                    p4est_topidx_t    which_tree,
                    int               num_outgoing,
                    p4est_quadrant_t *outgoing[],
                    int               num_incoming, 
                    p4est_quadrant_t *incoming[])
{
  QuadData_t          *parentData, *childData;

  /*--------------------------------------------------------
  | Coarsening -> Initialize new coarser quad from its
  |               children
  --------------------------------------------------------*/
  if (num_outgoing > 1)
  {
    int i, j, k;

    parentData = (QuadData_t *) incoming[0]->p.user_data;

    /*------------------------------------------------------
    | Init geometry and flow data for new quad
    ------------------------------------------------------*/
#ifdef P4_TO_P8
    init_quadGeomData3d(p4est, which_tree, 
                        incoming[0], parentData);
#else
    init_quadGeomData2d(p4est, which_tree, 
                        incoming[0], parentData);
#endif

    init_quadFlowData(parentData);

    /*------------------------------------------------------
    | Interpolate data from children
    ------------------------------------------------------*/
    for (i = 0; i < P4EST_CHILDREN; i++) 
    {
      childData = (QuadData_t *) outgoing[i]->p.user_data;

      for (j = 0; j < OCT_MAX_VARS; j++)
      {
        parentData->vars[j] += childData->vars[j];

        for (k = 0; k < P4EST_DIM; k++)
        {
          parentData->grad_vars[j][k] += childData->grad_vars[j][k];
        }
      }
    }

    /*------------------------------------------------------
    | Normalize with number of children
    ------------------------------------------------------*/
    for (j = 0; j < OCT_MAX_VARS; j++)
    {
      parentData->vars[j] /= P4EST_CHILDREN;

      for (k = 0; k < P4EST_DIM; k++)
      {
        parentData->grad_vars[j][k] /= P4EST_CHILDREN;
      }
    }

  }
  /*--------------------------------------------------------
  | Refinement -> Initialize new finer quads from their 
  |               parent
  --------------------------------------------------------*/
  else
  {
    parentData = (QuadData_t *) outgoing[0]->p.user_data;

    // Quad centroids
    octDouble *pxx = parentData->centroid;

    int i, j, k;

    for (i = 0; i < P4EST_CHILDREN; i++)
    {
      childData = (QuadData_t *) incoming[i]->p.user_data;

      /*----------------------------------------------------
      | Init geometry and flow data for new quad
      ----------------------------------------------------*/
#ifdef P4_TO_P8
      init_quadGeomData3d(p4est, which_tree, 
                          incoming[i], childData);
#else
      init_quadGeomData2d(p4est, which_tree, 
                          incoming[i], childData);
#endif

      init_quadFlowData(childData);

      /*----------------------------------------------------
      | Inpterpolate flow field data
      ----------------------------------------------------*/
      for (j = 0; j < OCT_MAX_VARS; j++)
      {
        childData->vars[j] = parentData->vars[j];

        octDouble *cxx = childData->centroid;

        for (k = 0; k < P4EST_DIM; k++)
        {
          const octDouble  dx = cxx[k] - pxx[k];

          childData->vars[j] += dx * parentData->grad_vars[j][k];

          childData->grad_vars[j][k] = parentData->grad_vars[j][k];
        }
      }
    }
  }

} /* interpQuadData() */

