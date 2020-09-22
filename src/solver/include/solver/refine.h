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
#ifndef SOLVER_REFINE_H
#define SOLVER_REFINE_H

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
* calcSqrErr_scalar()
*-----------------------------------------------------------
* This function estimates the approximation error 
* on a quadrant for the passive scalar variable.
*
* The estimate is computed by integrating the difference of a 
* constant approximation at the midpoint and a linear 
* approximation that interpolates at the midpoint.
*
* The square of the error estiamte for the state variables 
* contained is returned.
*
***********************************************************/
octDouble calcSqrErr(p4est_quadrant_t *q, int varIdx);

/***********************************************************
* refinement_scalarError()
*-----------------------------------------------------------
* Function to calculate the error estimate for the mesh
* refinement 
***********************************************************/
int refinement_scalarError(p4est_t          *p4est, 
                           p4est_topidx_t    which_tree,
                           p4est_quadrant_t *q);

/***********************************************************
* globalRefinement()
*-----------------------------------------------------------
* This is the general function to control the refinement
* of trees.
***********************************************************/
int globalRefinement(p4est_t          *p4est,
                     p4est_topidx_t    which_tree,
                     p4est_quadrant_t *q);


#endif /* SOLVER_REFINE_H */
