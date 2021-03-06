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
#ifndef SOLVER_SOLVER_H_INCLUDED
#define SOLVER_SOLVER_H_INCLUDED

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_vtk.h>
#include <p4est_iterate.h>
#else
#include <p8est_bits.h>
#include <p8est_vtk.h>
#endif

#include "solver/simData.h"

/***********************************************************
* solverRun()
*-----------------------------------------------------------
* Perform a transient simulation with initialized 
* simulation data
***********************************************************/
void solverRun(SimData_t *simData);


#endif /* SOLVER_SOLVER_H_INCLUDED */
