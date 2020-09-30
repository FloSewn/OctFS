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
#ifndef SOLVER_CONVECTIVEFLUX_H
#define SOLVER_CONVECTIVEFLUX_H

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#endif

#include "solver/typedefs.h"
#include "solver/simData.h"
#include "solver/quadData.h"

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
                      void *user_data);

#endif /* SOLVER_CONVECTIVEFLUX_H */
