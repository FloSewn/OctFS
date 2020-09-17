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
#ifndef SOLVER_SIMDATA_H
#define SOLVER_SIMDATA_H

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

#include "solver/typedefs.h"
#include "solver/util.h"


/***********************************************************
* Structure containing simulation parameters 
* (reference values, physical properties, constants...)
***********************************************************/
typedef struct SimParam_t
{
  /* Fluid viscosity */
  octDouble viscosity;

  /* Reference length */
  octDouble ref_length;
  /* Reference velocity */
  octDouble ref_velocity;
  /* Reference pressure */
  octDouble ref_pressure;

  /* User-defined initialization function for flow variables */
  octInitFun usrInitFun;
  /* User-defined refinement function */
  octRefineFun usrRefineFun;
  /* User-defined coarsening function */
  octCoarseFun usrCoarseFun;

} SimParam_t;

/***********************************************************
* Structure containing solver parameters 
* (max. number of iterations...)
***********************************************************/
typedef struct SolverParam_t
{
  octDouble epsilon;

  // Path to export directory
  char     *io_exportDir;
  // Prefix for export files
  char     *io_exportPrefix;

  // Number of quadrants per MPU
  octInt    nQuadMPU;
  // Minimum level of refinement for initialization
  octInt    minRefLvl;
  // Maximum level of refinement
  octInt    maxRefLvl;
  // Fill uniform
  octBool   fillUniform;
  // Recursive refinement
  octBool   recursive;
  // Re-Partition on coarsening
  octBool   partForCoarsen;

} SolverParam_t;

/***********************************************************
* Structure containing MPI parameters 
***********************************************************/
typedef struct MPIParam_t
{
  /* MPI Communicator */
  sc_MPI_Comm     mpiComm;

} MPIParam_t;

/***********************************************************
* Structure containing properties for the entire simulation
*   > Accessed through p4est->user_pointer
***********************************************************/
typedef struct SimData_t
{
  /* Simulation parameters */
  SimParam_t              *simParam;

  /* Solver parameters */
  SolverParam_t           *solverParam;

  /* MPI data */
  MPIParam_t              *mpiParam;

  /* p4est tree structure entry */
  p4est_t                 *p4est;

  /* p4est mesh connectivity */
  p4est_connectivity_t    *conn;

} SimData_t;

/***********************************************************
* init_simData()
*-----------------------------------------------------------
* Initializes the simulation data structure
***********************************************************/
SimData_t *init_simData(int          argc, 
                        char        *argv[], 
                        octInitFun   usrInitFun,
                        octRefineFun usrRefineFun,
                        octCoarseFun usrCoarseFun);

/***********************************************************
* init_simParam()
*-----------------------------------------------------------
* Initializes the simulation parameter structure
***********************************************************/
SimParam_t *init_simParam(octInitFun   usrInitFun,
                          octRefineFun usrRefineFun,
                          octCoarseFun usrCoarseFun);

/***********************************************************
* init_solverParam()
*-----------------------------------------------------------
* Initializes the solver parameter structure
***********************************************************/
SolverParam_t *init_solverParam(void);

/***********************************************************
* init_mpiParam()
*-----------------------------------------------------------
* Initializes the solver parameter structure
***********************************************************/
MPIParam_t *init_mpiParam(int argc, char *argv[]);


/***********************************************************
* destroy_simData()
*-----------------------------------------------------------
* Frees all memory of a SimData structure
***********************************************************/
void destroy_simData(SimData_t *simData);

/***********************************************************
* destroy_simParam()
*-----------------------------------------------------------
* Frees all memory of a SimParam structure
***********************************************************/
void destroy_simParam(SimParam_t *simParam);

/***********************************************************
* destroy_solverParam()
*-----------------------------------------------------------
* Frees all memory of a SolverParam structure
***********************************************************/
void destroy_solverParam(SolverParam_t *solverParam);

/***********************************************************
* destroy_mpiParam()
*-----------------------------------------------------------
* Frees all memory of a MPIParam structure
***********************************************************/
void destroy_mpiParam(MPIParam_t *mpiParam);

#endif /* SOLVER_SIMDATA_H */
