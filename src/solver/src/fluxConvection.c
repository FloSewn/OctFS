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
#include "solver/fluxConvection.h"
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
* Funtion to determine the upwind element. 
* var_out is the element with the outward facing normal
* var_in is the element with the inward facing normal.
***********************************************************/
#define UPWIND_DIR(mflux, var_out, var_in) \
  ( (mflux) > 0.0 ? (var_out) : (var_in) )



/***********************************************************
* addFlux_conv_imp()
*-----------------------------------------------------------
* Function to add the implicit part of convective fluxes.
* 
*   -> p4est_iter_face_t callback function
*
*-----------------------------------------------------------
* Arguments:
* *info      : Information about this quadrant that has 
*              been populated by the p4est_iterate()
* *user_data : user_data that is given to p4est_iterate,
*              in this case, it points to the ghost_data
*              array, which contains the SimVars_t data
*              for all of the ghost cells, that has been 
*              populated by p4est_ghost_exchange_data
*
***********************************************************/
void addFlux_conv_imp(p4est_iter_face_info_t *info,
                      void *user_data)
{
  p4est_t    *p4est    = info->p4est;
  SimData_t  *simData  = (SimData_t *) p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  octDouble fluxFac = simParam->tmp_fluxFac;
  int       varIdx  = simParam->tmp_varIdx;

  QuadData_t      *qData;
  QuadData_t      *ghostData = (QuadData_t *) user_data;

  sc_array_t *sides = &(info->sides);
  P4EST_ASSERT(sides->elem_count == 2);

  int i;

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
      qData = &ghostData[side[1]->is.full.quadid];
    else
      qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    const octDouble var_1 = qData->vars[varIdx];
    octDouble      *Ax_1  = qData->Ax_p;


    /*-----------------------------------------------------
    | Side 0: -> Take massfluxes
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      if (side[0]->is.hanging.is_ghost[i])
        qData = &ghostData[side[0]->is.hanging.quadid[i]];
      else
        qData = (QuadData_t*) side[0]->is.hanging.quad[i]->p.user_data;

      const octDouble mflux = qData->mflux[iface];
      const octDouble var_0 = qData->vars[varIdx];
      octDouble      *Ax_0  = qData->Ax_p;

      /*---------------------------------------------------
      | Determine upwind direction
      | -> var_0 has the outward facing normal
      |--------------------------------------------------*/
      const octDouble var_u = UPWIND_DIR(mflux,var_0,var_1);

      /*---------------------------------------------------
      | Add fluxes
      |--------------------------------------------------*/
      const octDouble flux = fluxFac * var_u * mflux;

      Ax_0[varIdx] += flux;
      Ax_1[varIdx] -= flux;

    } /* for (i = 0; i < P4EST_HALF; i++) */
  }
  else if (!side[0]->is_hanging && side[1]->is_hanging)
  {
    const int iface = side[1]->face;

    /*-----------------------------------------------------
    | Side 0: Large face 
    |----------------------------------------------------*/
    if (side[0]->is.full.is_ghost)
      qData = &ghostData[side[0]->is.full.quadid];
    else
      qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    const octDouble var_0 = qData->vars[varIdx];
    octDouble      *Ax_0  = qData->Ax_p;

    /*-----------------------------------------------------
    | Side 1: -> Take massfluxes
    | Hanging face: There are 2^(d-1) (P4EST_HALF) subfaces
    |----------------------------------------------------*/
    for (i = 0; i < P4EST_HALF; i++)
    {
      if (side[1]->is.hanging.is_ghost[i])
        qData = &ghostData[side[1]->is.hanging.quadid[i]];
      else
        qData = (QuadData_t*) side[1]->is.hanging.quad[i]->p.user_data;

      const octDouble mflux = qData->mflux[iface];
      const octDouble var_1 = qData->vars[varIdx];
      octDouble      *Ax_1  = qData->Ax_p;

      /*---------------------------------------------------
      | Determine upwind direction
      | -> var_1 has the outward facing normal
      |--------------------------------------------------*/
      const octDouble var_u = UPWIND_DIR(mflux,var_1,var_0);

      /*---------------------------------------------------
      | Add fluxes
      |--------------------------------------------------*/
      const octDouble flux = fluxFac * var_u * mflux;

      Ax_0[varIdx] -= flux;
      Ax_1[varIdx] += flux;

    } /* for (i = 0; i < P4EST_HALF; i++) */
  }
  else if (!side[0]->is_hanging && !side[1]->is_hanging)
  {
    const int iface = side[0]->face;

    /*-----------------------------------------------------
    | Side 0: Large face 
    |----------------------------------------------------*/
    if (side[0]->is.full.is_ghost)
      qData = &ghostData[side[0]->is.full.quadid];
    else
      qData = (QuadData_t *) side[0]->is.full.quad->p.user_data;

    const octDouble mflux = qData->mflux[iface];
    const octDouble var_0 = qData->vars[varIdx];
    octDouble      *Ax_0  = qData->Ax_p;

    /*-----------------------------------------------------
    | Side 1: Large face 
    |----------------------------------------------------*/
    if (side[1]->is.full.is_ghost)
      qData = &ghostData[side[1]->is.full.quadid];
    else
      qData = (QuadData_t *) side[1]->is.full.quad->p.user_data;

    const octDouble var_1 = qData->vars[varIdx];
    octDouble      *Ax_1  = qData->Ax_p;

    /*-----------------------------------------------------
    | Determine upwind direction
    | -> var_0 has the outward facing normal
    |----------------------------------------------------*/
    const octDouble var_u = UPWIND_DIR(mflux,var_0,var_1);

    /*-----------------------------------------------------
    | Add fluxes
    |----------------------------------------------------*/
    const octDouble flux = fluxFac * var_u * mflux;

    Ax_0[varIdx] += flux;
    Ax_1[varIdx] -= flux;
  }
  else
  {
    octPrint("UNDEFINED BEHAVIOUR!!!\n");
  }

} /* addFlux_conv_imp() */



