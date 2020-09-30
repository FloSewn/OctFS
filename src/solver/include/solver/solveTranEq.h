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
#ifndef SOLVER_SOLVETRANEQ_H
#define SOLVER_SOLVETRANEQ_H

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_iterate.h>
#endif

#include "solver/simData.h"

/***********************************************************
* resetSolverBuffers()
*-----------------------------------------------------------
* Function sets all solver buffers Ax, b, res for every 
* quad to zero
***********************************************************/
void resetSolverBuffers(p4est_iter_volume_info_t *info,
                        void *user_data);

/***********************************************************
* compute_b_tranEq()
*-----------------------------------------------------------
* This function sums up the right hand side of the equation
* system
*   Ax = b 
* that underlies a discretized transport equation.
***********************************************************/
void compute_b_tranEq(SimData_t *simData, int varIdx);

/***********************************************************
* solveTranEq()
*-----------------------------------------------------------
* Function to solve a transport equation for a specified
* variable <varIdx>.
***********************************************************/
void solveTranEq(SimData_t *simData, int varIdx);

#endif /* SOLVER_SOLVETRANEQ_H */
