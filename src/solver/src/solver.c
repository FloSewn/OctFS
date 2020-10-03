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
#include "solver/solver.h"
#include "solver/simData.h"
#include "solver/quadData.h"
#include "solver/dataIO.h"
#include "solver/refine.h"
#include "solver/coarsen.h"
#include "solver/gradients.h"
#include "solver/projection.h"
#include "solver/massflux.h"

#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#endif


/***********************************************************
* solverRun()
*-----------------------------------------------------------
* Perform a transient simulation with initialized 
* simulation data
***********************************************************/
void solverRun(SimData_t *simData)
{
  SimParam_t *simParam        = simData->simParam;
  SolverParam_t *solverParam  = simData->solverParam;

  int step;
  octDouble time;

  int refinePeriod      = solverParam->refinePeriod;
  int repartitionPeriod = solverParam->repartitionPeriod;
  int writePeriod       = solverParam->writePeriod;

  octBool adaptGrid     = solverParam->adaptGrid;

  octDouble dt         = simParam->timestep;
  octDouble simTimeTot = simParam->simTimeTot;

  /*--------------------------------------------------------
  | Initialize gradients
  --------------------------------------------------------*/
  int idx;
  for (idx = 0; idx < OCT_MAX_VARS; idx++)
    computeGradients(simData, idx); 

  /*--------------------------------------------------------
  | The main loop
  --------------------------------------------------------*/
  for (time = 0.0, step = 0; 
       time < simTimeTot; time += dt, step += 1)
  {
    simParam->simTime += dt;

    octPrint("TIME STEP %d", step);

    /*------------------------------------------------------
    | Perform a refinement of the domain
    ------------------------------------------------------*/
    if (  !(step % refinePeriod) 
        && (step > 0) 
        && (adaptGrid == TRUE) )
    {
      p4est_refine_ext(simData->p4est,
                       solverParam->recursive,
                       solverParam->maxRefLvl,
                       globalRefinement,
                       NULL,
                       interpQuadData);

      p4est_coarsen_ext(simData->p4est, 
                        solverParam->recursive, 
                        0,
                        globalCoarsening, 
                        NULL,
                        interpQuadData);

      p4est_balance_ext(simData->p4est, 
                        P4EST_CONNECT_FACE, 
                        NULL,
                        interpQuadData);

      p4est_ghost_destroy(simData->ghost);
      P4EST_FREE(simData->ghostData);
      simData->ghost = NULL;
      simData->ghostData = NULL;
    }

    /*------------------------------------------------------
    | Repartition domain
    |-----------------------------------------------------*/
    if (    (step > 0)
        && !(step % repartitionPeriod) 
        &&  (adaptGrid == TRUE) ) 
    {
      p4est_partition(simData->p4est, 
                      solverParam->partForCoarsen, 
                      NULL);

      if (simData->ghost) 
      {
        p4est_ghost_destroy(simData->ghost);
        P4EST_FREE(simData->ghostData);
        simData->ghost = NULL;
        simData->ghostData = NULL;
      }
    }

    /*------------------------------------------------------
    | Synchronize ghost data
    |-----------------------------------------------------*/
    if (!simData->ghost) {
      simData->ghost = p4est_ghost_new(simData->p4est, 
                                       P4EST_CONNECT_FULL);
      simData->ghostData = P4EST_ALLOC(QuadData_t, 
                          simData->ghost->ghosts.elem_count);
      p4est_ghost_exchange_data(simData->p4est, 
                                simData->ghost, 
                                simData->ghostData);
    }

    /*------------------------------------------------------
    | Solve projection step
    |-----------------------------------------------------*/
    doProjectionStep(simData);

    /*------------------------------------------------------
    | Print out solver data to user
    |-----------------------------------------------------*/

    /*------------------------------------------------------
    | Write solution
    |-----------------------------------------------------*/
    if (!(step % writePeriod)) 
    {
      octPrint("WRITE SOLUTION FILE FOR STEP %d", step);
      writeSolutionVtk(simData, step);
    }


  } /* for (time, step ...) */

  /*------------------------------------------------------
  | Write final solution
  |-----------------------------------------------------*/
  octPrint("WRITE SOLUTION FILE FOR STEP %d", step);
  writeSolutionVtk(simData, step);

  /*--------------------------------------------------------
  | Release ghost data
  |-------------------------------------------------------*/
  P4EST_FREE(simData->ghostData);
  p4est_ghost_destroy(simData->ghost);


} /* solverRun() */
