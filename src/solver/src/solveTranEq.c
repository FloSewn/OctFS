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
* resetSolverBuffers_Ax()
*-----------------------------------------------------------
* Function sets all solver buffers Ax for every 
* quad to zero and lets Ax_p point to Ax.
*
*   -> p4est_iter_cell_t callback function
***********************************************************/
void resetSolverBuffers_Ax(p4est_iter_volume_info_t *info,
                           void *user_data)
{
  p4est_quadrant_t   *q = info->quad;

  QuadData_t     *quadData = (QuadData_t *) q->p.user_data;

  int i; 
  for (i = 0; i < OCT_MAX_VARS; i++)
    quadData->Ax[i]  = 0.0;

  quadData->Ax_p = quadData->Ax;

} /* resetSolverBuffers_Ax() */

/***********************************************************
* resetSolverBuffers_b()
*-----------------------------------------------------------
* Function sets all solver buffers b for every 
* quad to zero and lets Ax_p point to b.
*
*   -> p4est_iter_cell_t callback function
***********************************************************/
void resetSolverBuffers_b(p4est_iter_volume_info_t *info,
                           void *user_data)
{
  p4est_quadrant_t   *q = info->quad;

  QuadData_t     *quadData = (QuadData_t *) q->p.user_data;

  int i; 
  for (i = 0; i < OCT_MAX_VARS; i++)
    quadData->b[i]   = 0.0;

  quadData->Ax_p = quadData->b;

} /* resetSolverBuffers_b() */

/***********************************************************
* compute_b_tranEq()
*-----------------------------------------------------------
* This function sums up the right hand side of the equation
* system
*   Ax = b 
* that underlies a discretized transport equation.
***********************************************************/
void compute_b_tranEq(SimData_t *simData, int varIdx)
{
  /*--------------------------------------------------------
  | Compute flux factors for chosen discretization scheme
  --------------------------------------------------------*/
  SimParam_t *simParam = simData->simParam;

  int scheme            = simParam->tempScheme;
  simParam->tmp_fluxFac = simParam->tempFluxFac[scheme]-1.0;
  simParam->tmp_varIdx  = varIdx;

  /*--------------------------------------------------------
  | Update gradient
  --------------------------------------------------------*/
  computeGradients(simData, varIdx); 

  /*--------------------------------------------------------
  | Empty buffer and set pointer Ax_p to b
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                resetSolverBuffers_b,  // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback*/

  /*--------------------------------------------------------
  | Add convective fluxes
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est,      
                simData->ghost,      
                (void *) simData->ghostData, 
                NULL,             // quad callback
                addFlux_conv_imp, // face callback
#ifdef P4_TO_P8
                NULL,             // edge callback
#endif
                NULL);            // corner callback

  /*--------------------------------------------------------
  | Add diffusive fluxes
  --------------------------------------------------------*/

  /*--------------------------------------------------------
  | Add sources
  --------------------------------------------------------*/

  /*--------------------------------------------------------
  | Add temporal derivative terms
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                addTimeDerivative,     // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback

} /* compute_b_tranEq() */

/***********************************************************
* compute_Ax_tranEq()
*-----------------------------------------------------------
* This function sums up the left hand side of the equation
* system
*   Ax = b 
* that underlies a discretized transport equation.
***********************************************************/
void compute_Ax_tranEq(SimData_t *simData, int varIdx)
{
  /*--------------------------------------------------------
  | Compute flux factors for chosen discretization scheme
  --------------------------------------------------------*/
  SimParam_t *simParam = simData->simParam;

  int scheme            = simParam->tempScheme;
  simParam->tmp_fluxFac = simParam->tempFluxFac[scheme];
  simParam->tmp_varIdx  = varIdx;

  /*--------------------------------------------------------
  | Update gradient
  --------------------------------------------------------*/
  computeGradients(simData, varIdx); 

  /*--------------------------------------------------------
  | Empty buffer and set pointer Ax_p to b
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                resetSolverBuffers_Ax, // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback*/

  /*--------------------------------------------------------
  | Add convective fluxes
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est,      
                simData->ghost,      
                (void *) simData->ghostData, 
                NULL,             // quad callback
                addFlux_conv_imp, // face callback
#ifdef P4_TO_P8
                NULL,             // edge callback
#endif
                NULL);            // corner callback

  /*--------------------------------------------------------
  | Add diffusive fluxes
  --------------------------------------------------------*/

  /*--------------------------------------------------------
  | Add sources
  --------------------------------------------------------*/

  /*--------------------------------------------------------
  | Add temporal derivative terms
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                addTimeDerivative,     // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback

} /* compute_b_tranEq() */

/***********************************************************
* solveTranEq()
*-----------------------------------------------------------
* Function to solve a transport equation for a specified
* variable <varIdx>.
***********************************************************/
void solveTranEq(SimData_t *simData, int varIdx)
{
  /*--------------------------------------------------------
  | Compute right hand side b
  --------------------------------------------------------*/
  compute_b_tranEq(simData, varIdx);

  /*--------------------------------------------------------
  | Solve transport equation using Krylov solver
  --------------------------------------------------------*/

} /* solveTranEq() */
