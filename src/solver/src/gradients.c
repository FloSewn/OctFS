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
#include "solver/quadData.h"
#include "solver/simData.h"
#include "solver/util.h"
#include "aux/dbg.h"

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#endif


/***********************************************************
* Constants  
***********************************************************/
static int varIdx;

/***********************************************************
* resetDerivatives()
*-----------------------------------------------------------
* Function to reset the derivatives of a quadrant 
*   -> p4est_iter_volume_t callback function
***********************************************************/
void resetDerivatives(p4est_iter_volume_info_t *info,
                      void *user_data)
{
  p4est_quadrant_t   *q = info->quad;

  QuadData_t     *quadData = (QuadData_t *) q->p.user_data;
  QuadFlowData_t *flowData = &(quadData->flowData);

  int i;

  for (i = 0; i < P4EST_DIM; i++)
  {
    flowData->grad_vars[varIdx][i] = 0.0;
  }

} /* resetDerivatives() */

/***********************************************************
* divideByVolume()
*-----------------------------------------------------------
* Function to divide by volume for Green-Gauss gradients.
*   -> p4est_iter_volume_t callback function
***********************************************************/
void divideByVolume(p4est_iter_volume_info_t *info,
                    void                     *user_data)
{
  p4est_quadrant_t   *q = info->quad;

  QuadData_t     *quadData = (QuadData_t *) q->p.user_data;
  QuadFlowData_t *flowData = &(quadData->flowData);
  QuadGeomData_t *geomData = &(quadData->geomData);

  int i;

  octDouble vol = 1.0 / geomData->volume;

  for (i = 0; i < P4EST_DIM; i++)
  {
    flowData->grad_vars[varIdx][i] *= vol;
  }

} /* divideByVolume() */

/***********************************************************
* computeGradGauss()
*-----------------------------------------------------------
* Function to calculate the spatial gradient 
* based on a Green-Gauss algrithm.
*   -> p4est_iter_face_t callback function
***********************************************************/
void computeGradGauss(p4est_iter_face_info_t *info,
                      void                   *user_data)
{
  int i, j;

  sc_array_t *sides = &(info->sides);
  P4EST_ASSERT(sides->elem_count == 2);

  QuadData_t      *ghostData = (QuadData_t *) user_data;
  QuadData_t      *qData;
  QuadFlowData_t  *fData;
  QuadGeomData_t  *gData;

  p4est_quadrant_t *quad;

  octDouble  var_0, var_1;
  octDouble  *grad_0, *grad_1;

  /*-------------------------------------------------------
  | No boundaries implemented yet -> every face has two 
  | sides
  |------------------------------------------------------*/
  p4est_iter_face_side_t *side[2];
  side[0] = p4est_iter_fside_array_index_int(sides, 0);
  side[1] = p4est_iter_fside_array_index_int(sides, 1);

  /*-------------------------------------------------------
  | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
  |------------------------------------------------------*/
  if (side[0]->is_hanging && !side[1]->is_hanging)
  {
    const int iface = side[0]->face;

    /*-----------------------------------------------------
    | Side 1: Normal face 
    |----------------------------------------------------*/
    quad = side[1]->is.full.quad;  

    if (side[1]->is.full.is_ghost)
      qData = &ghostData[side[1]->is.full.quadid];
    else
      qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    fData  = &(qData->flowData);
    var_1  = fData->vars[varIdx];
    grad_1 = fData->grad_vars[varIdx];

    /*-----------------------------------------------------
    | Side 0: -> Set normals
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      quad = side[0]->is.hanging.quad[i];

      if (side[0]->is.hanging.is_ghost[i])
        qData = &ghostData[side[0]->is.hanging.quadid[i]];
      else
        qData = (QuadData_t*) side[0]->is.hanging.quad[i]->p.user_data;

      fData = &(qData->flowData);
      gData = &(qData->geomData);

      var_0  = fData->vars[varIdx];
      grad_0 = fData->grad_vars[varIdx];

      const octDouble nx = gData->normals[iface][0];
      const octDouble ny = gData->normals[iface][1];
#ifdef P4_TO_P8
      const octDouble nz = gData->normals[iface][2];
#endif

      /*---------------------------------------------------
      | Add flux contribution
      |--------------------------------------------------*/
      const octDouble var_f = 0.5 * (var_0 + var_1);

      const octDouble grad_x = nx * var_f;
      const octDouble grad_y = ny * var_f;
#ifdef P4_TO_P8
      const octDouble grad_z = nz * var_f;
#endif

      grad_0[0] += grad_x;
      grad_1[0] -= grad_x;

      grad_0[1] += grad_y;
      grad_1[1] -= grad_y;

#ifdef P4_TO_P8
      grad_0[2] += grad_z;
      grad_1[2] -= grad_z;
#endif
    }
  }
  else if (!side[0]->is_hanging && side[1]->is_hanging)
  {
    const int iface = side[1]->face;

    /*-----------------------------------------------------
    | Side 0:  Normal face 
    |----------------------------------------------------*/
    quad = side[0]->is.full.quad;  

    if (side[0]->is.full.is_ghost)
      qData = &ghostData[side[0]->is.full.quadid];
    else
      qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    fData = &(qData->flowData);

    var_0  = fData->vars[varIdx];
    grad_0 = fData->grad_vars[varIdx];

    /*-----------------------------------------------------
    | Side 1:
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    | -> Set normals
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      quad = side[1]->is.hanging.quad[i];

      if (side[1]->is.hanging.is_ghost[i])
        qData = &ghostData[side[1]->is.hanging.quadid[i]];
      else
        qData = (QuadData_t*) side[1]->is.hanging.quad[i]->p.user_data;

      fData = &(qData->flowData);
      gData = &(qData->geomData);

      var_1   = fData->vars[varIdx];
      grad_1  = fData->grad_vars[varIdx];

      const octDouble nx = gData->normals[iface][0];
      const octDouble ny = gData->normals[iface][1];
#ifdef P4_TO_P8
      const octDouble nz = gData->normals[iface][2];
#endif

      /*---------------------------------------------------
      | Add flux contribution
      |--------------------------------------------------*/
      const octDouble var_f = 0.5 * (var_0 + var_1);
      
      const octDouble grad_x = nx * var_f;
      const octDouble grad_y = ny * var_f;
#ifdef P4_TO_P8
      const octDouble grad_z = nz * var_f;
#endif

      grad_0[0] -= grad_x;
      grad_1[0] += grad_x;

      grad_0[1] -= grad_y;
      grad_1[1] += grad_y;

#ifdef P4_TO_P8
      grad_0[2] -= grad_z;
      grad_1[2] += grad_z;
#endif
    }
  }
  else if (!side[0]->is_hanging && !side[1]->is_hanging)
  {
    const int iface = side[0]->face;

    /*-----------------------------------------------------
    | Side 0:
    | Normal face 
    |----------------------------------------------------*/
    quad = side[0]->is.full.quad;  

    if (side[0]->is.full.is_ghost)
      qData = &ghostData[side[0]->is.full.quadid];
    else
      qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    fData = &(qData->flowData);
    gData = &(qData->geomData);

    var_0  = fData->vars[varIdx];
    grad_0 = fData->grad_vars[varIdx];

    const octDouble nx = gData->normals[iface][0];
    const octDouble ny = gData->normals[iface][1];
    const octDouble nz = gData->normals[iface][2];

    /*-----------------------------------------------------
    | Side 1:
    | Normal face 
    |----------------------------------------------------*/
    quad = side[1]->is.full.quad;  

    if (side[1]->is.full.is_ghost)
      qData = &ghostData[side[1]->is.full.quadid];
    else
      qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    fData = &(qData->flowData);

    var_1  = fData->vars[varIdx];
    grad_1 = fData->grad_vars[varIdx];

    /*-----------------------------------------------------
    | Add flux contribution
    |----------------------------------------------------*/
    const octDouble var_f = 0.5 * (var_0 + var_1);

    const octDouble grad_x = nx * var_f;
    const octDouble grad_y = ny * var_f;
#ifdef P4_TO_P8
    const octDouble grad_z = nz * var_f;
#endif

    grad_0[0] += grad_x;
    grad_1[0] -= grad_x;

    grad_0[1] += grad_y;
    grad_1[1] -= grad_y;

#ifdef P4_TO_P8
    grad_0[2] += grad_z;
    grad_1[2] -= grad_z;
#endif

  }
  else
  {
    printf("UNDEFINED BEHAVIOUR!!!\n");
  }

} /* computeGradGauss() */


/***********************************************************
* computeGradient()
*-----------------------------------------------------------
* Function to calculate the spatial gradient within the 
* domain.
***********************************************************/
void computeGradients(p4est_t       *p4est, 
                      p4est_ghost_t *ghost,
                      QuadData_t    *ghost_data,
                      int            varIdx_)
{
  /*-------------------------------------------------------
  | Set variable index for gradient calculation
  -------------------------------------------------------*/
  varIdx = varIdx_;

  /*-------------------------------------------------------
  | Green-Gauss gradient estimation
  -------------------------------------------------------*/
  p4est_iterate(p4est, 
                ghost, 
                (void *) ghost_data,
                resetDerivatives,  // cell callback
                computeGradGauss,   // face callback
#ifdef P4_TO_P8
                NULL,     // edge callback
#endif
                NULL);   // corner callback*/

  /*-------------------------------------------------------
  | Scaling by volume
  -------------------------------------------------------*/
  p4est_iterate(p4est, 
                ghost, 
                (void *) ghost_data,
                divideByVolume,  // cell callback
                NULL,     // face callback
#ifdef P4_TO_P8
                NULL,     // edge callback
#endif
                NULL);   // corner callback


} /* computeGradients(...) */
