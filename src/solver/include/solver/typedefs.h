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
#ifndef SOLVER_TYPEDEFS_H
#define SOLVER_TYPEDEFS_H

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
* double and int length 
***********************************************************/
#define octDouble double
#define octInt    int
#define octBool   int

/***********************************************************
* Solver variables
***********************************************************/
#define OCT_MAX_VARS       17 /* Max. number of variables */
#define OCT_SOLVER_VARS    11 /* Number of solver variab. */
#define OCT_VARNAME_LENGTH 32 /* Max. var. name length    */

#define QUAD_BUF_VARS    10 /* No. of lin. solver buffs.*/
#define PARAM_BUF_VARS   10 /* No. of lin. solver buffs.*/

/***********************************************************
* Solver indices
*-----------------------------------------------------------
* Indices for buffer variables of linear solver are
* stored for every quad structure prior to field variables
***********************************************************/
typedef enum 
{
  SAX,  /* Holds results for the product Ax               */
  SB,   /* Holds results for the right hand side b        */
  SVN,  /* Holds new updated values of variable data      */
  SR,   /* Contains solution to (b - Ax) at iteration n   */
  SR0,  /* Contains solution to (b - Ax) at iteration 0   */
  SP,   /* Direction for new solution                     */
  SV,   /* Holds A*p                                      */
  SH,   /* x - alpha * p                                  */
  SS,   /* r - alpha * v                                  */
  ST,   /*                                                */
  SRES, /* buffer for general (b - Ax)                    */
  IRHO,
  IVX,
  IVY,
  IVZ,
  IP,
  IS
} VarIndex;

/***********************************************************
* Temporal schemes
***********************************************************/
typedef enum 
{ 
  EULER_EXPLICIT,
  EULER_IMPLICIT,
  CRANK_NICOLSON
} TempScheme;

/***********************************************************
* Indices for buffer variables of linear solver that are
* stored in the simParam struct
***********************************************************/
typedef enum
{
  PR0,  /* rho_0                                          */
  PA,   /* alpha                                          */
  PO,   /* omega                                          */
  PR,   /* rho                                            */
  PB,   /* beta                                           */
  PRES, /* residual                                       */
  PGRES /* global residual                                */
} SimParamBufIndex;

/***********************************************************
* Typedefs for simData.h
***********************************************************/
typedef struct SimParam_t       SimParam_t;
typedef struct SolverParam_t    SolverParam_t;
typedef struct MPIParam_t       MPIParam_t;
typedef struct SimData_t        SimData_t;

/***********************************************************
* Typedefs for quadData.h
***********************************************************/
typedef struct QuadData_t       QuadData_t;

/***********************************************************
* Initialization function pointer for user 
***********************************************************/
typedef void (*octInitFun)   (QuadData_t *quadData);
typedef int  (*octRefineFun) (p4est_t          *p4est,
                              p4est_topidx_t    which_tree,
                              p4est_quadrant_t *q);
typedef int  (*octCoarseFun) (p4est_t          *p4est,
                              p4est_topidx_t    which_tree,
                              p4est_quadrant_t *children[]);


#endif /* SOLVER_TYPEDEFS_H */
