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
* linSolve_calcGlobResidual()
*-----------------------------------------------------------
* Linear solver function for the calculaiton of the global 
* residual.
* The local residual is always stored in vars[SRES] and
* the global residual is stored in sbuf[PRES]
***********************************************************/
octDouble linSolve_calcGlobResidual(SimData_t *simData,
                                    computeAx  cmpAx,
                                    int         xId,
                                    int        AxId,
                                    int         bId)
{
  SimParam_t *simParam = simData->simParam;

  const int n_elements  = simData->p4est->global_num_quadrants;
  const octDouble n_inv = 1. / (octDouble) n_elements;

  /*--------------------------------------------------------
  | vars[AxId] = A * vars[xId]
  | reset simParam->tmp_xId to xId, since its changed 
  | inside cmpAx()
  --------------------------------------------------------*/
  int xId_old = simParam->tmp_xId;
  cmpAx(simData, xId, AxId);
  simParam->tmp_xId = xId_old;

  /*--------------------------------------------------------
  | vars[SRES] = (1.0)*vars[bId] + (-1.0)*vars[AxId] 
  --------------------------------------------------------*/
  linSolve_fieldSum(simData, bId, AxId, SRES, 1.0,  -1.0);

  /*--------------------------------------------------------
  | sbuf[PRES] = sum( vars[SRES] * vars[SRES] )
  --------------------------------------------------------*/
  linSolve_scalarProd(simData, SRES, SRES, PRES);
  simParam->sbuf[PRES] = n_inv * sqrt(simParam->sbuf[PRES]);

  return simParam->sbuf[PRES];

} /* linSolve_calcGlobResidual() */

/***********************************************************
* linSolve_exchangeScalarBuffer()
*-----------------------------------------------------------
* Linear solver function to sum a scalar buffer 
* variable over all MPI processes 
***********************************************************/
void linSolve_exchangeScalarBuffer(SimData_t      *simData, 
                                   int             sbufId,
                                   sc_MPI_Datatype varType,
                                   sc_MPI_Op       mpiType)
{
  SimParam_t *simParam = simData->simParam;

  octDouble buf = simParam->sbuf[sbufId];

  sc_MPI_Allreduce(&buf, &simParam->sbuf[sbufId],
                   1, varType, mpiType,
                   simData->mpiParam->mpiComm);

} /* linSolve_exchangeScalarBuffer() */

/***********************************************************
* linSolve_scalarProd()
*-----------------------------------------------------------
* Linear solver function to multiply two field variables 
* a and b according to
*
*   c = a_i * b_i
*
* and store the in a scalar variable c.
***********************************************************/
void linSolve_scalarProd(SimData_t *simData, 
                         int aId, int bId, int cId)
{
  int prodBuf[3] = { aId, bId, cId };

  SimParam_t *simParam = simData->simParam;

  simParam->sbuf[cId] = 0.0;

  p4est_iterate(simData->p4est, NULL, (void *) prodBuf,
                linSolve_scalarProd_cb, // cell callback
                NULL,                   // face callback
#ifdef P4_TO_P8
                NULL,                   // edge callback
#endif
                NULL);                  // corner callback*/

  /*--------------------------------------------------------
  | Exchange data among all processes
  --------------------------------------------------------*/
  linSolve_exchangeScalarBuffer(simData, cId, 
                                  sc_MPI_DOUBLE, 
                                  sc_MPI_SUM);

} /* linSolve_scalarProd() */

/***********************************************************
* linSolve_fieldProd()
*-----------------------------------------------------------
* Linear solver function to multiply two field variables 
* a and b according to
*
*   c_i = a_i * b_i
*
* and store the in another field variable c.
***********************************************************/
void linSolve_fieldProd(SimData_t *simData, 
                        int aId, int bId, int cId)
{
  int prodBuf[3] = { aId, bId, cId };

  p4est_iterate(simData->p4est, NULL, (void *) prodBuf,
                linSolve_fieldProd_cb, // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback*/

} /* linSolve_fieldProd() */

/***********************************************************
* linSolve_scalarSum()
*-----------------------------------------------------------
* Linear solver function to add two field variables 
* a and b according to
*
*   c = w_a * a_i + w_b * b_i
*
* and store the in a scalar variable c.
***********************************************************/
void linSolve_scalarSum(SimData_t *simData, 
                        int aId, int bId, int cId,
                        octDouble w_a, octDouble w_b)
{
  octDouble sumBuf[5] = { (octDouble) aId,
                          (octDouble) bId,
                          (octDouble) cId,
                           w_a, w_b };

  SimParam_t *simParam = simData->simParam;

  simParam->sbuf[cId] = 0.0;

  p4est_iterate(simData->p4est, NULL, (void *) sumBuf,
                linSolve_scalarSum_cb, // cell callback
                NULL,                 // face callback
#ifdef P4_TO_P8
                NULL,                 // edge callback
#endif
                NULL);                // corner callback*/

} /* linSolve_scalarSum() */

/***********************************************************
* linSolve_fieldSum()
*-----------------------------------------------------------
* Linear solver function to add two field variables a and b
* according to
*
*   c_i = w_a * a_i + w_b * b_i
*
* and store the in another field variable c.
***********************************************************/
void linSolve_fieldSum(SimData_t *simData, 
                       int aId, int bId, int cId, 
                       octDouble w_a, octDouble w_b)
{
  octDouble sumBuf[5] = { (octDouble) aId,
                          (octDouble) bId,
                          (octDouble) cId,
                           w_a, w_b };

  p4est_iterate(simData->p4est, NULL, (void *) sumBuf,
                linSolve_fieldSum_cb, // cell callback
                NULL,                 // face callback
#ifdef P4_TO_P8
                NULL,                 // edge callback
#endif
                NULL);                // corner callback*/

} /* linSolve_fieldSum() */

/***********************************************************
* linSolve_fieldCopy()
*-----------------------------------------------------------
* Linear solver function to add two field variables a into
* a field variable b
***********************************************************/
void linSolve_fieldCopy(SimData_t *simData, 
                       int aId, int bId)
{
  int cpyBuf[2] = {  aId, bId };

  p4est_iterate(simData->p4est, NULL, (void *) cpyBuf,
                linSolve_fieldCopy_cb, // cell callback
                NULL,                  // face callback
#ifdef P4_TO_P8
                NULL,                  // edge callback
#endif
                NULL);                 // corner callback*/

} /* linSolve_fieldCopy() */



/***********************************************************
* linSolve_fieldSum_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_fieldSum()
***********************************************************/
void linSolve_fieldSum_cb(p4est_iter_volume_info_t *info,
                          void *user_data)
{
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;

  int aId = (int)((octDouble *) user_data)[0];
  int bId = (int)((octDouble *) user_data)[1];
  int cId = (int)((octDouble *) user_data)[2];

  octDouble w_a = ((octDouble *) user_data)[3];
  octDouble w_b = ((octDouble *) user_data)[4];

  octDouble a = quadData->vars[aId];
  octDouble b = quadData->vars[bId];

  quadData->vars[cId] = w_a*a + w_b*b;

} /* linSolve_fieldSum_cb() */

/***********************************************************
* linSolve_scalarSum_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_scalarSum()
***********************************************************/
void linSolve_scalarSum_cb(p4est_iter_volume_info_t *info,
                           void *user_data)
{
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;
  SimData_t  *simData  = (SimData_t*)info->p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  int aId = (int)((octDouble *) user_data)[0];
  int bId = (int)((octDouble *) user_data)[1];
  int cId = (int)((octDouble *) user_data)[2];

  octDouble w_a = ((octDouble *) user_data)[3];
  octDouble w_b = ((octDouble *) user_data)[4];

  octDouble a = quadData->vars[aId];
  octDouble b = quadData->vars[bId];

  simParam->sbuf[cId] += w_a*a + w_b*b;

} /* linSolve_scalarSum_cb() */

/***********************************************************
* linSolve_fieldProd_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_fieldProd()
***********************************************************/
void linSolve_fieldProd_cb(p4est_iter_volume_info_t *info,
                           void *user_data)
{
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;

  int aId = ((int *) user_data)[0];
  int bId = ((int *) user_data)[1];
  int cId = ((int *) user_data)[2];

  octDouble a = quadData->vars[aId];
  octDouble b = quadData->vars[bId];

  quadData->vars[cId] = a * b;

} /* linSolve_fieldProd_cb() */

/***********************************************************
* linSolve_scalarProd_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_scalarProd()
***********************************************************/
void linSolve_scalarProd_cb(p4est_iter_volume_info_t *info,
                            void *user_data)
{
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;
  SimData_t  *simData  = (SimData_t*)info->p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  int aId = ((int *) user_data)[0];
  int bId = ((int *) user_data)[1];
  int cId = ((int *) user_data)[2];

  octDouble a = quadData->vars[aId];
  octDouble b = quadData->vars[bId];

  simParam->sbuf[cId] += a * b;

} /* linSolve_scalarProd_cb() */

/***********************************************************
* linSolve_fieldCopy_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_fieldCopy()
***********************************************************/
void linSolve_fieldCopy_cb(p4est_iter_volume_info_t *info,
                           void *user_data)
{
  QuadData_t *quadData = (QuadData_t*)info->quad->p.user_data;

  int aId = ((int *) user_data)[0];
  int bId = ((int *) user_data)[1];

  quadData->vars[bId] = quadData->vars[aId];

} /* linSolve_fieldCopy_cb() */



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

  SimParam_t *simParam = simData->simParam;
  int         xId      = simParam->tmp_xId;

  const octDouble vol     = quadData->volume;
  const octDouble dt      = simParam->timestep;
  const octDouble rho     = quadData->vars[IRHO];
  const octDouble b       = quadData->vars[SB];

  quadData->vars[xId] = b * dt / vol / rho;

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
                       int        xId)
{
  SimParam_t *simParam  = simData->simParam;
  int n_elements        = simData->p4est->global_num_quadrants;
  const octDouble n_inv = 1. / (octDouble) n_elements;

  int k = 0;

  /*--------------------------------------------------------
  | Init scalar solver buffers
  --------------------------------------------------------*/
  simParam->sbuf[PR0]   = 1.0;
  simParam->sbuf[PA]    = 1.0;
  simParam->sbuf[PO]    = 1.0;
  simParam->sbuf[PR]    = 0.0;
  simParam->sbuf[PB]    = 0.0;
  simParam->sbuf[PRES]  = 0.0;
  simParam->sbuf[PGRES] = 0.0;

  /*--------------------------------------------------------
  | Threshold parameters
  --------------------------------------------------------*/
  int kMin = 2;
  int kMax = 50;

  octDouble eps = 1e-6;

  /*--------------------------------------------------------
  | Compute new Ax
  --------------------------------------------------------*/
  cmpAx(simData, xId, SAX);

  /*--------------------------------------------------------
  | vars[SR] = (1.0)*vars[SB] + (-1.0)*vars[SAX]
  --------------------------------------------------------*/
  linSolve_fieldSum(simData, SB, SAX, SR, 1.0, -1.0);

  /*--------------------------------------------------------
  | vars[SR0] = vars[SR] 
  --------------------------------------------------------*/
  linSolve_fieldCopy(simData, SR, SR0);

  /*--------------------------------------------------------
  | sbuf[PGRES] = sqrt( sum( vars[SR] * vars[SR] ) ) / N
  --------------------------------------------------------*/
  linSolve_scalarProd(simData, SR, SR, PGRES);
  simParam->sbuf[PGRES] = n_inv * sqrt(simParam->sbuf[PGRES]);

  while( k < kMax )
  {
    k++;

    /*------------------------------------------------------
    | sbuf[PR] = sum( vars[SR0] * vars[SR] )
    ------------------------------------------------------*/
    linSolve_scalarProd(simData, SR0, SR, PR);

    /*------------------------------------------------------
    | Update buffers beta (PB) and rho_0 (PR0) 
    ------------------------------------------------------*/
    octDouble rho   = simParam->sbuf[PR];
    octDouble rho_0 = simParam->sbuf[PR0];
    octDouble alpha = simParam->sbuf[PA];
    octDouble omega = simParam->sbuf[PO];

    simParam->sbuf[PB]  = (rho / (SMALL+rho_0)) 
                        * (alpha / (SMALL+omega));
    simParam->sbuf[PR0] = rho;

    /*------------------------------------------------------
    | 1) vars[SP] = (1.0)*vars[SP] + (-sbuf[PO])*vars[SV]
    | 2) vars[SP] = (1.0)*vars[SR] + ( sbuf[PB])*vars[SP]
    ------------------------------------------------------*/
    linSolve_fieldSum(simData, SP, SV, SP, 
                      1.0, -simParam->sbuf[PO]);
    linSolve_fieldSum(simData, SR, SP, SP, 
                      1.0,  simParam->sbuf[PB]);

    /*------------------------------------------------------
    | Compute v = A*p
    | reset simParam->tmp_xId to xId, since its changed 
    | inside cmpAx()
    ------------------------------------------------------*/
    cmpAx(simData, SP, SV);
    simParam->tmp_xId = xId;

    /*------------------------------------------------------
    | sbuf[PA] = sum( vars[SR0] * vars[SV] )
    ------------------------------------------------------*/
    linSolve_scalarProd(simData, SR0, SV, PA);

    alpha = 1. / (SMALL + simParam->sbuf[PA]);
    simParam->sbuf[PA] = alpha * simParam->sbuf[PR]; 

    /*------------------------------------------------------
    | vars[SH] = (1.0)*vars[xId] + (sbuf[PA])*vars[SP]
    ------------------------------------------------------*/
    linSolve_fieldSum(simData, xId, SP, SH, 
                      1.0,  simParam->sbuf[PA]);

    /*------------------------------------------------------
    | Calculate global residual for the equation system
    | A*h = b and store it in sbuf[PRES]
    ------------------------------------------------------*/
    linSolve_calcGlobResidual(simData, cmpAx, SH, SAX, SB);

    /*------------------------------------------------------
    | Check if vars[SH] is accuarte enough
    | if yes -> set as new solution and resume
    ------------------------------------------------------*/
    if ( simParam->sbuf[PRES] < eps && k > kMin )
    {
      linSolve_fieldCopy(simData, SH, xId);
      break;
    }

    /*------------------------------------------------------
    | vars[SS] = (1.0)*vars[SR] + (-sbuf[PA])*vars[SV]
    ------------------------------------------------------*/
    linSolve_fieldSum(simData, SR, SV, SS, 
                      1.0,  -simParam->sbuf[PA]);

    /*------------------------------------------------------
    | vars[ST] = A*vars[SS]
    ------------------------------------------------------*/
    cmpAx(simData, SS, ST);
    simParam->tmp_xId = xId;

    /*------------------------------------------------------
    | sbuf[PO] = sum( vars[ST] * vars[ST] )
    ------------------------------------------------------*/
    linSolve_scalarProd(simData, ST, ST, PO);
    omega = 1. / (simParam->sbuf[PO] + SMALL);

    /*------------------------------------------------------
    | sbuf[PO] = sum( vars[ST] * vars[SS] )
    ------------------------------------------------------*/
    linSolve_scalarProd(simData, ST, SS, PO);
    simParam->sbuf[PO] *= omega;

    /*------------------------------------------------------
    | vars[xId] = (1.0)*vars[SH] + (sbuf[PO])*vars[SS]
    ------------------------------------------------------*/
    linSolve_fieldSum(simData, SH, SS, xId, 
                      1.0,  simParam->sbuf[PO]);

    /*------------------------------------------------------
    | Calculate global residual for the equation system
    | A*h = b and store it in sbuf[PRES]
    ------------------------------------------------------*/
    linSolve_calcGlobResidual(simData, cmpAx, xId, SAX, SB);

    /*------------------------------------------------------
    | Check if vars[xId] is accuarte enough
    ------------------------------------------------------*/
    if ( simParam->sbuf[PRES] < eps && k > kMin )
      break;

    /*------------------------------------------------------
    | vars[SR] = (1.0)*vars[SS] + (-sbuf[PO])*vars[ST]
    ------------------------------------------------------*/
    linSolve_fieldSum(simData, SS, ST, SR, 
                      1.0, -simParam->sbuf[PO]);

  } /* while( k < kMax ) */


} /* linSolve_bicgstab() */


/***********************************************************
* solve_explicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an explicit method.
***********************************************************/
void solve_explicit_sequential(SimData_t *simData, 
                               int        xId)
{
  SimParam_t *simParam = simData->simParam;
  simParam->tmp_xId    = xId;

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
                               int        xId)
{
  SimParam_t *simParam = simData->simParam;
  simParam->tmp_xId    = xId;

  /*--------------------------------------------------------
  | Solve linear equation system using Krylov solver
  --------------------------------------------------------*/
  linSolve_bicgstab(simData, cmpAx, xId);


} /* solve_implicit_sequential()*/
