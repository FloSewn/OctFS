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
#ifndef SOLVER_LINEARSOLVER_H
#define SOLVER_LINEARSOLVER_H

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
* computeAx()
*-----------------------------------------------------------
* Function pointer to function for calculation of Ax.
***********************************************************/
typedef void (*computeAx) (SimData_t *simData,
                           int        varIdx,
                           int        sbufIdx);

/***********************************************************
* resetSolverBuffers()
*-----------------------------------------------------------
* Function sets all solver buffer variables to zero.
*
*   -> p4est_iter_volume_t callback function
***********************************************************/
void resetSolverBuffers(p4est_iter_volume_info_t *info,
                        void *user_data);

/***********************************************************
* addRightHandSide()
*-----------------------------------------------------------
* Function adds the right hand side b to the
* solution 
*
*   -> p4est_iter_cell_t callback function
***********************************************************/
void addRightHandSide(p4est_iter_volume_info_t *info,
                      void *user_data);

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
                       int        varIdx);

/***********************************************************
* solve_explicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an explicit method.
***********************************************************/
void solve_explicit_sequential(SimData_t *simData, 
                               int        varIdx);

/***********************************************************
* solve_implicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an implicit method.
***********************************************************/
void solve_implicit_sequential(SimData_t *simData, 
                               computeAx  cmpAx,
                               int        varIdx);

#endif /* SOLVER_LINEARSOLVER_H */
