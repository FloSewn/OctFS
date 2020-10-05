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
static int gradVarIdx = -1;


/***********************************************************
* setGradVarIdx()
*-----------------------------------------------------------
* Function to set the index of the gradient variable to
* compute.
***********************************************************/
static void setGradVarIdx(int varIdx)
{
  gradVarIdx = varIdx;
}

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

  int i;

  for (i = 0; i < P4EST_DIM; i++)
  {
    quadData->grad_vars[gradVarIdx][i] = 0.0;
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

  int i;

  octDouble vol = 1.0 / quadData->volume;

  for (i = 0; i < P4EST_DIM; i++)
  {
    quadData->grad_vars[gradVarIdx][i] *= vol;
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
  int i;

  QuadData_t      *qDat0, *qDat1;
  QuadData_t      *ghostData = (QuadData_t *) user_data;

  sc_array_t *sides = &(info->sides);
  P4EST_ASSERT(sides->elem_count == 2);

  /*-------------------------------------------------------
  | every face has two sides
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
    | Side 1: Large face 
    |----------------------------------------------------*/
    if (side[1]->is.full.is_ghost)
      qDat1 = &ghostData[side[1]->is.full.quadid];
    else
      qDat1 = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    /*-----------------------------------------------------
    | Side 0: -> Set normals
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      if (side[0]->is.hanging.is_ghost[i])
        qDat0 = &ghostData[side[0]->is.hanging.quadid[i]];
      else
        qDat0 = (QuadData_t*) side[0]->is.hanging.quad[i]->p.user_data;


      /*---------------------------------------------------
      | Add flux contribution
      |--------------------------------------------------*/
      const octDouble nx = qDat0->normals[iface][0];
      const octDouble ny = qDat0->normals[iface][1];
#ifdef P4_TO_P8
      const octDouble nz = qDat0->normals[iface][2];
#endif
      const octDouble var_1  = qDat1->vars[gradVarIdx];
      const octDouble var_0  = qDat0->vars[gradVarIdx];
      const octDouble var_f = 0.5 * (var_0 + var_1);

      const octDouble grad_x = nx * var_f;
      const octDouble grad_y = ny * var_f;
#ifdef P4_TO_P8
      const octDouble grad_z = nz * var_f;
#endif

      qDat0->grad_vars[gradVarIdx][0] += grad_x;
      qDat1->grad_vars[gradVarIdx][0] -= grad_x;

      qDat0->grad_vars[gradVarIdx][1] += grad_y;
      qDat1->grad_vars[gradVarIdx][1] -= grad_y;

#ifdef P4_TO_P8
      qDat0->grad_vars[gradVarIdx][2] += grad_z;
      qDat1->grad_vars[gradVarIdx][2] -= grad_z;
#endif
    }
  }
  else if (!side[0]->is_hanging && side[1]->is_hanging)
  {
    const int iface = side[1]->face;

    /*-----------------------------------------------------
    | Side 0:  Large face 
    |----------------------------------------------------*/
    if (side[0]->is.full.is_ghost)
      qDat0 = &ghostData[side[0]->is.full.quadid];
    else
      qDat0 = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    /*-----------------------------------------------------
    | Side 1:
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    | -> Set normals
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      if (side[1]->is.hanging.is_ghost[i])
        qDat1 = &ghostData[side[1]->is.hanging.quadid[i]];
      else
        qDat1 = (QuadData_t*) side[1]->is.hanging.quad[i]->p.user_data;

      /*---------------------------------------------------
      | Add flux contribution
      |--------------------------------------------------*/
      const octDouble nx = qDat1->normals[iface][0];
      const octDouble ny = qDat1->normals[iface][1];
#ifdef P4_TO_P8
      const octDouble nz = qDat1->normals[iface][2];
#endif
      const octDouble var_1  = qDat1->vars[gradVarIdx];
      const octDouble var_0  = qDat0->vars[gradVarIdx];
      const octDouble var_f = 0.5 * (var_0 + var_1);
      
      const octDouble grad_x = nx * var_f;
      const octDouble grad_y = ny * var_f;
#ifdef P4_TO_P8
      const octDouble grad_z = nz * var_f;
#endif

      qDat0->grad_vars[gradVarIdx][0] -= grad_x;
      qDat1->grad_vars[gradVarIdx][0] += grad_x;

      qDat0->grad_vars[gradVarIdx][1] -= grad_y;
      qDat1->grad_vars[gradVarIdx][1] += grad_y;

#ifdef P4_TO_P8
      qDat0->grad_vars[gradVarIdx][2] -= grad_z;
      qDat1->grad_vars[gradVarIdx][2] += grad_z;
#endif
    }
  }
  else if (!side[0]->is_hanging && !side[1]->is_hanging)
  {
    const int iface = side[0]->face;

    /*-----------------------------------------------------
    | Side 0: Large face 
    |----------------------------------------------------*/
    if (side[0]->is.full.is_ghost)
      qDat0 = &ghostData[side[0]->is.full.quadid];
    else
      qDat0 = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    /*-----------------------------------------------------
    | Side 1:
    | Normal face 
    |----------------------------------------------------*/
    if (side[1]->is.full.is_ghost)
      qDat1 = &ghostData[side[1]->is.full.quadid];
    else
      qDat1 = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    /*-----------------------------------------------------
    | Add flux contribution
    |----------------------------------------------------*/
    const octDouble nx = qDat0->normals[iface][0];
    const octDouble ny = qDat0->normals[iface][1];
#ifdef P4_TO_P8
    const octDouble nz = qDat0->normals[iface][2];
#endif
    const octDouble var_0  = qDat0->vars[gradVarIdx];
    const octDouble var_1  = qDat1->vars[gradVarIdx];
    const octDouble var_f = 0.5 * (var_0 + var_1);

    const octDouble grad_x = nx * var_f;
    const octDouble grad_y = ny * var_f;
#ifdef P4_TO_P8
    const octDouble grad_z = nz * var_f;
#endif

    qDat0->grad_vars[gradVarIdx][0] += grad_x;
    qDat1->grad_vars[gradVarIdx][0] -= grad_x;

    qDat0->grad_vars[gradVarIdx][1] += grad_y;
    qDat1->grad_vars[gradVarIdx][1] -= grad_y;

#ifdef P4_TO_P8
    qDat0->grad_vars[gradVarIdx][2] += grad_z;
    qDat1->grad_vars[gradVarIdx][2] -= grad_z;
#endif

  }
  else
  {
    octPrint("UNDEFINED BEHAVIOUR!!!\n");
  }

} /* computeGradGauss() */


/***********************************************************
* computeGradient()
*-----------------------------------------------------------
* Function to calculate the spatial gradient within the 
* domain.
***********************************************************/
void computeGradients(SimData_t *simData, int varIdx)
{
  p4est_t       *p4est       = simData->p4est;
  p4est_ghost_t *ghost       = simData->ghost;
  QuadData_t    *ghostData   = simData->ghostData;
  
  /*-------------------------------------------------------
  | Set variable index for gradient calculation
  -------------------------------------------------------*/
  setGradVarIdx(varIdx);

  /*--------------------------------------------------------
  | Exchange data
  --------------------------------------------------------*/
  p4est_ghost_exchange_data(simData->p4est, 
                            simData->ghost, 
                            simData->ghostData);

  /*-------------------------------------------------------
  | Green-Gauss gradient estimation
  -------------------------------------------------------*/
  p4est_iterate(p4est, 
                ghost, 
                (void *) ghostData,
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
                (void *) ghostData,
                divideByVolume,  // cell callback
                NULL,     // face callback
#ifdef P4_TO_P8
                NULL,     // edge callback
#endif
                NULL);   // corner callback



} /* computeGradients(...) */
