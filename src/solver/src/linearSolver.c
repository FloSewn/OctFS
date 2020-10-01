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
#include "solver/solveTranEq.h"
#include "solver/simData.h"
#include "solver/quadData.h"
#include "solver/gradients.h"
#include "solver/fluxConvection.h"
#include "solver/timeIntegral.h"

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
* addRightHandSide()
*-----------------------------------------------------------
* Function adds the right hand side b to the
* solution 
*
*   -> p4est_iter_cell_t callback function
***********************************************************/
void addRightHandSide(p4est_iter_volume_info_t *info,
                      void *user_data)
{
  p4est_quadrant_t  *q = info->quad;

  p4est_t    *p4est    = info->p4est;
  SimData_t  *simData  = (SimData_t *) p4est->user_pointer;
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;

  SimParam_t *simParam  = simData->simParam;
  int         varIdx    = simParam->tmp_varIdx;


  octDouble vol     = quadData->volume;
  octDouble dt      = simParam->timestep;
  octDouble rho     = quadData->vars[IRHO];
  octDouble fac     = dt / (vol * rho);

  quadData->vars[varIdx] = quadData->b[varIdx] * fac;

} /* resetSolverBuffers_Ax() */


/***********************************************************
* solve_explicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an explicit method.
***********************************************************/
void solve_explicit_sequential(SimData_t *simData, 
                               int        varIdx)
{
  SimParam_t *simParam  = simData->simParam;
  simParam->tmp_varIdx  = varIdx;

  /*--------------------------------------------------------
  | Add right hand side to solution 
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                addRightHandSide,  // cell callback
                NULL,              // face callback
#ifdef P4_TO_P8
                NULL,              // edge callback
#endif
                NULL);             // corner callback*/

} /* solve_explicit_sequential() */
