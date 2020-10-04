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
#include "solver/timeIntegral.h"
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
* addTimeDerivative()
*-----------------------------------------------------------
* Function to add the temporal derivative.
* 
*   -> p4est_iter_volume_t callback function
***********************************************************/
void addTimeDerivative(p4est_iter_volume_info_t *info,
                       void *user_data)
{
  SimData_t  *simData  = (SimData_t*)info->p4est->user_pointer;
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;
  SimParam_t *simParam = simData->simParam;

  int       xId  = simParam->tmp_xId;
  int       AxId = simParam->tmp_AxId;
  octDouble vol  = quadData->volume;
  octDouble dt   = simParam->timestep;
  octDouble rho  = quadData->vars[IRHO];
  octDouble var  = quadData->vars[xId];

  quadData->vars[AxId] += vol * var * rho / dt; 

} /* addTimeDerivative() */
