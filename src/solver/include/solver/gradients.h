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
#ifndef SOLVER_GRADIENTS_H
#define SOLVER_GRADIENTS_H

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
#include "solver/quadData.h"

/***********************************************************
* setGradVarIdx()
*-----------------------------------------------------------
* Function to set the index of the gradient variable to
* compute.
***********************************************************/
static void setGradVarIdx(int varIdx);

/***********************************************************
* resetDerivatives()
*-----------------------------------------------------------
* Function to reset the derivatives of a quadrant 
*   -> p4est_iter_volume_t callback function
***********************************************************/
void resetDerivatives(p4est_iter_volume_t *info,
                      void *user_data);

/***********************************************************
* divideByVolume()
*-----------------------------------------------------------
* Function to divide by volume for Green-Gauss gradients.
*   -> p4est_iter_volume_t callback function
***********************************************************/
void divideByVolume(p4est_iter_volume_info_t *info,
                    void                     *user_data);

/***********************************************************
* computeGradGauss()
*-----------------------------------------------------------
* Function to calculate the spatial gradient 
* based on a Green-Gauss algrithm.
*   -> p4est_iter_face_t callback function
***********************************************************/
void computeGradGauss(p4est_iter_face_info_t *info,
                      void                   *user_data);

/***********************************************************
* computeGradient()
*-----------------------------------------------------------
* Function to calculate the spatial gradient within the 
* domain.
***********************************************************/
void computeGradients(SimData_t *simData, int varIdx);

#endif /* SOLVER_GRADIENTS_H */
