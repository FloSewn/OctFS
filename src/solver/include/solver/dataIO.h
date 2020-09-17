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
#ifndef SOLVER_DATAIO_H
#define SOLVER_DATAIO_H

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
#include "solver/simData.h"


/***********************************************************
* Index for writing data
***********************************************************/
static int io_idx;

/***********************************************************
* Function to concatenate two strings.
* Appends <str_2> to the end of <str_1>.
*-----------------------------------------------------------
* Arguments:
* *str_1, str_2 : strings to concatenate
*
***********************************************************/
char *concat_string(const char *str_1, const char *str_2);

/***********************************************************
* Callback function for interpolating the solution from 
* quadrant midpoints to corners.
*
* The function p4est_iterate() takes as an argument a 
* p4est_iter_volume_t callback function, which it executes at 
* every local quadrant (see p4est_iterate.h).  
* This function matches the p4est_iter_volume_t prototype.
*
* Use the callback function to interpolate the state
* variable to the corners, and write those corners into an 
* array so that they can be written out.
*-----------------------------------------------------------
* Arguments:
* *info       : the information about this quadrant that 
*               has been  populated by p4est_iterate()
* user_data   : the user_data that was given as an argument to
*               p4est_iterate: in this case, it points to the
*               array of corner values that we want to write.
*               The values for the corner of the quadrant
*               described by a info are written during the
*               execution of the callback.
*
***********************************************************/
static void interpolate_solution(p4est_iter_volume_info_t *info, 
                                 void *user_data);

/***********************************************************
* writeSolutionVtk()
*-----------------------------------------------------------
* Function to write the solution of a single timestep 
* to vtk format
***********************************************************/
void writeSolutionVtk(SimData_t *simData, 
                      int        step, 
                      int        varIdx);
                      

#endif /* SOLVER_DATAIO_H */
