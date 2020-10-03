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
#define OCT_MAX_VARS        6 /* Max. number of variables */
#define OCT_VARNAME_LENGTH 32 /* Max. var. name length    */

#define SOLVER_BUF_VARS    10 /* No. of lin. solver buffs.*/

/***********************************************************
* Solver indices
***********************************************************/
typedef enum 
{
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
* Indices for buffer variables of linear solver 
***********************************************************/
typedef enum 
{
  LS_AX,  /* Holds results for the product Ax             */
  LS_B,   /* Holds results for the right hand side b      */
  LS_VN,  /* Holds new updated values of variable data    */
  LS_R,   /* Contains solution to (b - Ax) at iteration n */
  LS_R0,  /* Contains solution to (b - Ax) at iteration 0 */
  LS_P,   /* Direction for new solution                   */
  LS_V,   /* Holds A*p                                    */
  LS_H,   /* x - alpha * p                                */
  LS_S,   /* r - alpha * v                                */
  LS_t    /*                                              */
} LinSolveIndex;

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
