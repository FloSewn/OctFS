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
#include "solver/massflux.h"
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
* resetMasflux()
*-----------------------------------------------------------
* Function sets all mass fluxes at element interfaces to 
* zero
*
*   -> p4est_iter_volume_t callback function
***********************************************************/
void resetMassflux(p4est_iter_face_info_t *info,
                   void                   *user_data)
{
  QuadData_t      *qData;
  QuadData_t      *ghostData = (QuadData_t *) user_data;

  sc_array_t *sides = &(info->sides);
  P4EST_ASSERT(sides->elem_count == 2);

  p4est_iter_face_side_t *side[2];
  side[0] = p4est_iter_fside_array_index_int(sides, 0);
  side[1] = p4est_iter_fside_array_index_int(sides, 1);

  const int iface_0 = side[0]->face;
  const int iface_1 = side[1]->face;

  /*--------------------------------------------------------
  | Side 0
  |-------------------------------------------------------*/
  if (side[0]->is.full.is_ghost)
    qData = &ghostData[side[0]->is.full.quadid];
  else
    qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

  qData->mflux[iface_0] = 0.0;

  /*--------------------------------------------------------
  | Side 1
  |-------------------------------------------------------*/
  if (side[1]->is.full.is_ghost)
    qData = &ghostData[side[1]->is.full.quadid];
  else
    qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

  qData->mflux[iface_1] = 0.0;

} /* resetSolverBuffers() */

/***********************************************************
* computeMassflux()
*-----------------------------------------------------------
* Function to reconstruct the mass fluxes at element 
* interfaces. The massfluxes are always stored at the 
* smaller face, if two elements of different quad-size
* are adjacent. For similar quad-sizes, the lower 
* indexed element is chosen.
*
*   -> p4est_iter_face_t callback function
***********************************************************/
void computeMassflux(p4est_iter_face_info_t *info,
                     void                   *user_data)
{
  QuadData_t      *qData;
  QuadData_t      *ghostData = (QuadData_t *) user_data;

  sc_array_t *sides = &(info->sides);
  P4EST_ASSERT(sides->elem_count == 2);


  int i;

  /*-------------------------------------------------------
  | No boundaries implemented yet -> every face has two 
  | sides
  |------------------------------------------------------*/
  p4est_iter_face_side_t *side[2];
  side[0] = p4est_iter_fside_array_index_int(sides, 0);
  side[1] = p4est_iter_fside_array_index_int(sides, 1);

  const int iface_0 = side[0]->face;
  const int iface_1 = side[1]->face;

  /*-------------------------------------------------------
  | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
  |------------------------------------------------------*/
  if (side[0]->is_hanging && !side[1]->is_hanging)
  {
    /*-----------------------------------------------------
    | Side 1: Normal face 
    |----------------------------------------------------*/
    if (side[1]->is.full.is_ghost)
      qData = &ghostData[side[1]->is.full.quadid];
    else
      qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    const octDouble u1 = qData->vars[IVX];
    const octDouble v1 = qData->vars[IVY];
#ifdef P4_TO_P8
    const octDouble w1 = qData->vars[IVZ];
#endif

    octDouble mflux_1 = qData->mflux[iface_1];

    /*-----------------------------------------------------
    | Side 0: -> Set normals
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      if (side[0]->is.hanging.is_ghost[i])
        qData = &ghostData[side[0]->is.hanging.quadid[i]];
      else
        qData = (QuadData_t*) side[0]->is.hanging.quad[i]->p.user_data;

      const octDouble u0 = qData->vars[IVX];
      const octDouble v0 = qData->vars[IVY];
#ifdef P4_TO_P8
      const octDouble w0 = qData->vars[IVZ];
#endif

      const octDouble nx = qData->normals[iface_0][0];
      const octDouble ny = qData->normals[iface_0][1];
#ifdef P4_TO_P8
      const octDouble nz = qData->normals[iface_0][2];
#endif

      /*---------------------------------------------------
      | Compute mass flux
      |--------------------------------------------------*/
      const octDouble um = 0.5 * (u0 + u1);
      const octDouble vm = 0.5 * (v0 + v1);
#ifdef P4_TO_P8
      const octDouble wm = 0.5 * (w0 + w1);
#endif

#ifdef P4_TO_P8
      const octDouble mflux = nx * um + ny * vm + nz * wm;
#else
      const octDouble mflux = nx * um + ny * vm; 
#endif

      /*---------------------------------------------------
      | Add mass fluxes to faces
      |--------------------------------------------------*/
      octDouble mflux_0 = qData->mflux[iface_0];

      mflux_0 += mflux;
      mflux_1 -= mflux;

    } /* for (i = 0; i < P4EST_HALF; i++) */
      
  } 
  else if (!side[0]->is_hanging && side[1]->is_hanging)
  {
    /*-----------------------------------------------------
    | Side 0: Normal face 
    |----------------------------------------------------*/
    if (side[0]->is.full.is_ghost)
      qData = &ghostData[side[0]->is.full.quadid];
    else
      qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    const octDouble u0 = qData->vars[IVX];
    const octDouble v0 = qData->vars[IVY];
#ifdef P4_TO_P8
    const octDouble w0 = qData->vars[IVZ];
#endif

    octDouble mflux_0 = qData->mflux[iface_0];

    /*-----------------------------------------------------
    | Side 1: -> Set normals
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      if (side[1]->is.hanging.is_ghost[i])
        qData = &ghostData[side[1]->is.hanging.quadid[i]];
      else
        qData = (QuadData_t*) side[1]->is.hanging.quad[i]->p.user_data;

      const octDouble u1 = qData->vars[IVX];
      const octDouble v1 = qData->vars[IVY];
#ifdef P4_TO_P8
      const octDouble w1 = qData->vars[IVZ];
#endif

      const octDouble nx = qData->normals[iface_1][0];
      const octDouble ny = qData->normals[iface_1][1];
#ifdef P4_TO_P8
      const octDouble nz = qData->normals[iface_1][2];
#endif

      /*---------------------------------------------------
      | Compute mass flux
      |--------------------------------------------------*/
      const octDouble um = 0.5 * (u0 + u1);
      const octDouble vm = 0.5 * (v0 + v1);
#ifdef P4_TO_P8
      const octDouble wm = 0.5 * (w0 + w1);
#endif

#ifdef P4_TO_P8
      const octDouble mflux = nx * um + ny * vm + nz * wm;
#else
      const octDouble mflux = nx * um + ny * vm; 
#endif

      /*---------------------------------------------------
      | Add mass fluxes to faces
      |--------------------------------------------------*/
      octDouble mflux_1 = qData->mflux[iface_1];

      mflux_0 -= mflux;
      mflux_1 += mflux;

    } /* for (i = 0; i < P4EST_HALF; i++) */
  }
  else if (!side[0]->is_hanging && !side[1]->is_hanging)
  {
    /*-----------------------------------------------------
    | Side 0: Normal face 
    |----------------------------------------------------*/
    if (side[0]->is.full.is_ghost)
      qData = &ghostData[side[0]->is.full.quadid];
    else
      qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    const octDouble u0 = qData->vars[IVX];
    const octDouble v0 = qData->vars[IVY];
#ifdef P4_TO_P8
    const octDouble w0 = qData->vars[IVZ];
#endif

    octDouble mflux_0 = qData->mflux[iface_0];

    /*-----------------------------------------------------
    | Use normals from side 0
    |----------------------------------------------------*/
    const octDouble nx = qData->normals[iface_0][0];
    const octDouble ny = qData->normals[iface_0][1];
#ifdef P4_TO_P8
    const octDouble nz = qData->normals[iface_0][2];
#endif


    /*-----------------------------------------------------
    | Side 1: Normal face 
    |----------------------------------------------------*/
    if (side[1]->is.full.is_ghost)
      qData = &ghostData[side[1]->is.full.quadid];
    else
      qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    const octDouble u1 = qData->vars[IVX];
    const octDouble v1 = qData->vars[IVY];
#ifdef P4_TO_P8
    const octDouble w1 = qData->vars[IVZ];
#endif

    octDouble mflux_1 = qData->mflux[iface_1];


    /*---------------------------------------------------
    | Compute mass flux
    |--------------------------------------------------*/
    const octDouble um = 0.5 * (u0 + u1);
    const octDouble vm = 0.5 * (v0 + v1);
#ifdef P4_TO_P8
    const octDouble wm = 0.5 * (w0 + w1);
#endif

#ifdef P4_TO_P8
    const octDouble mflux = nx * um + ny * vm + nz * wm;
#else
    const octDouble mflux = nx * um + ny * vm; 
#endif

    /*---------------------------------------------------
    | Add mass fluxes to faces
    |--------------------------------------------------*/
    mflux_0 += mflux;
    mflux_1 -= mflux;

  } 
  else
  {
    octPrint("UNDEFINED BEHAVIOUR!!!\n");
  }

} /* computeMassflux() */

/***********************************************************
* initMassfluxes()
*-----------------------------------------------------------
* Compute the massfluxes at element boundaries.
***********************************************************/
void initMassfluxes(SimData_t *simData)
{
  p4est_t       *p4est       = simData->p4est;
  p4est_ghost_t *ghost       = simData->ghost;
  QuadData_t    *ghostData   = simData->ghostData;

  /*-------------------------------------------------------
  | Set all massfluxes to zero
  -------------------------------------------------------*/
  p4est_iterate(p4est, 
                ghost, 
                (void *) ghostData,
                NULL,           // cell callback
                resetMassflux,  // face callback
#ifdef P4_TO_P8
                NULL,           // edge callback
#endif
                NULL);          // corner callback*/

  /*-------------------------------------------------------
  | Reconstruct mass fluxes at interfaces
  -------------------------------------------------------*/
  p4est_iterate(p4est, 
                ghost, 
                (void *) ghostData,
                NULL,             // cell callback
                computeMassflux,  // face callback
#ifdef P4_TO_P8
                NULL,             // edge callback
#endif
                NULL);            // corner callback*/

} /* calcMassfluxes() */

