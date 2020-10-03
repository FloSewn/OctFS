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
#include "solver/linearSolver.h"

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
* calcResiualDirection()
*-----------------------------------------------------------
* Calculate the direction (b - Ax) and compute initial
* squared residual norm r^2 = sum_i{ (b_i-Ax_i)^2 }
*
*   -> p4est_iter_volume_t callback function
***********************************************************/
void calcResidualDirection(p4est_iter_volume_info_t *info,
                           void *user_data)
{
  p4est_quadrant_t  *q = info->quad;
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;

  p4est_t    *p4est    = info->p4est;
  SimData_t  *simData  = (SimData_t *) p4est->user_pointer;

  SimParam_t *simParam  = simData->simParam;
  int         varIdx    = simParam->tmp_varIdx;

  octDouble Ax = quadData->sbuf[LS_AX][varIdx];
  octDouble b  = quadData->sbuf[LS_B][varIdx];

  octDouble dir = b - Ax;

  quadData->sbuf[LS_R][varIdx]  = dir;
  quadData->sbuf[LS_R0][varIdx] = dir;

  simParam->tmp_globRes += dir * dir;

} /* calcResidualDirection() */


/***********************************************************
* resetSolverBuffers()
*-----------------------------------------------------------
* Function sets all solver buffer variables to zero.
*
*   -> p4est_iter_volume_t callback function
***********************************************************/
void resetSolverBuffers(p4est_iter_volume_info_t *info,
                        void *user_data)
{
  p4est_quadrant_t  *q = info->quad;
  QuadData_t *quadData = (QuadData_t *) q->p.user_data;

  int i,j;


  for (j = 0; j < OCT_MAX_VARS; j++)
  {
    for (i = 0; i < SOLVER_BUF_VARS; i++)
      quadData->sbuf[i][j] = 0.0;

    quadData->res[j] = 0.0;
  }


} /* resetSolverBuffers_Ax() */

/***********************************************************
* addRightHandSide()
*-----------------------------------------------------------
* Function adds the right hand side b to the
* solution 
*
*   -> p4est_iter_volume_t callback function
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

  quadData->vars[varIdx] = quadData->sbuf[LS_B][varIdx] * fac;

} /* addRightHandSide() */


/***********************************************************
* linSolve_bicgstab()
*-----------------------------------------------------------
* Iterative solver for an equation system 
*
*   A x = b
*
* using a biconjugate gradient stabilized method (BICGSTAB)
*
***********************************************************/
void linSolve_bicgstab(SimData_t *simData,
                       computeAx  cmpAx,
                       int        varIdx)
{
  SimParam_t *simParam  = simData->simParam;

  /*--------------------------------------------------------
  | Init scalars
  --------------------------------------------------------*/
  int i, k = 0;

  octDouble rho_0   = 1.0;
  octDouble alpha   = 1.0;
  octDouble omega   = 1.0;
  octDouble rho     = 1.0;
  octDouble beta    = 0.0;
  octDouble abs_res = 0.0;

  /*--------------------------------------------------------
  | Threshold parameters
  --------------------------------------------------------*/
  int kMin = 2;
  int kMax = 10;

  /*--------------------------------------------------------
  | Compute new Ax
  --------------------------------------------------------*/
  cmpAx(simData, varIdx, LS_AX);

  /*--------------------------------------------------------
  | Compute direction (b - Ax)
  | and compute initial residual norm
  --------------------------------------------------------*/
  simParam->tmp_globRes = 0.0;

  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                calcResidualDirection, // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback*/

  octPrint("GLOBAL RESIDUAL: %lf", simParam->tmp_globRes);



} /* linSolve_bicgstab() */


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


/***********************************************************
* solve_implicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an implicit method.
***********************************************************/
void solve_implicit_sequential(SimData_t *simData, 
                               computeAx  cmpAx,
                               int        varIdx)
{
  SimParam_t *simParam  = simData->simParam;
  simParam->tmp_varIdx  = varIdx;

  /*--------------------------------------------------------
  | Initialize Krylov solver buffer variables
  --------------------------------------------------------*/
  p4est_iterate(simData->p4est, simData->ghost, 
                (void *) simData->ghostData,
                resetSolverBuffers, // cell callback
                NULL,              // face callback
#ifdef P4_TO_P8
                NULL,              // edge callback
#endif
                NULL);             // corner callback*/

  /*--------------------------------------------------------
  | Solve linear equation system using Krylov solver
  --------------------------------------------------------*/
  linSolve_bicgstab(simData, cmpAx, varIdx);


} /* solve_implicit_sequential()*/
