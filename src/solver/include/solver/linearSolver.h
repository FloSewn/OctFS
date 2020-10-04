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
                           int        xId,
                           int        sbufIdx);

/***********************************************************
* linSolve_exchangeScalarBuffer()
*-----------------------------------------------------------
* Linear solver function to sum a scalar buffer 
* variable over all MPI processes 
***********************************************************/
void linSolve_exchangeScalarBuffer(SimData_t      *simData, 
                                   int             sbufId,
                                   sc_MPI_Datatype varType,
                                   sc_MPI_Op       mpiType);


/***********************************************************
* linSolve_scalarProd()
*-----------------------------------------------------------
* Linear solver function to multiply two field variables 
* a and b according to
*
*   c = a_i * b_i
*
* and store the in a scalar variable c.
***********************************************************/
void linSolve_scalarProd(SimData_t *simData, 
                         int aId, int bId, int cId);

/***********************************************************
* linSolve_fieldProd()
*-----------------------------------------------------------
* Linear solver function to multiply two field variables 
* a and b according to
* *   c_i = a_i * b_i
*
* and store the in another field variable c.
***********************************************************/
void linSolve_fieldProd(SimData_t *simData, 
                        int aId, int bId, int cId);


/***********************************************************
* linSolve_scalarSum()
*-----------------------------------------------------------
* Linear solver function to add two field variables 
* a and b according to
*
*   c = w_a * a_i + w_b * b_i
*
* and store the in a scalar variable c.
***********************************************************/
void linSolve_scalarSum(SimData_t *simData, 
                        int aId, int bId, int cId,
                        octDouble w_a, octDouble w_b);

/***********************************************************
* linSolve_fieldSum()
*-----------------------------------------------------------
* Linear solver function to add two field variables a and b
* according to
*
*   c_i = w_a * a_i + w_b * b_i
*
* and store the in another field variable c.
***********************************************************/
void linSolve_fieldSum(SimData_t *simData, 
                       int aId, int bId, int cId, 
                       octDouble w_a, octDouble w_b);

/***********************************************************
* linSolve_fieldCopy()
*-----------------------------------------------------------
* Linear solver function to add two field variables a into
* a field variable b
***********************************************************/
void linSolve_fieldCopy(SimData_t *simData, 
                       int aId, int bId);


/***********************************************************
* linSolve_fieldSum_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_fieldSum()
***********************************************************/
void linSolve_fieldSum_cb(p4est_iter_volume_info_t *info,
                          void *user_data);

/***********************************************************
* linSolve_scalarSum_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_scalarSum()
***********************************************************/
void linSolve_scalarSum_cb(p4est_iter_volume_info_t *info,
                           void *user_data);

/***********************************************************
* linSolve_fieldProd_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_fieldProd()
***********************************************************/
void linSolve_fieldProd_cb(p4est_iter_volume_info_t *info,
                           void *user_data);

/***********************************************************
* linSolve_scalarProd_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_scalarProd()
***********************************************************/
void linSolve_scalarProd_cb(p4est_iter_volume_info_t *info,
                            void *user_data);

/***********************************************************
* linSolve_fieldCopy_cb()
*-----------------------------------------------------------
* p4est_iter_volume_t callback function for 
* linSolve_fieldCopy()
***********************************************************/
void linSolve_fieldCopy_cb(p4est_iter_volume_info_t *info,
                           void *user_data);







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
                       int        xId);

/***********************************************************
* solve_explicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an explicit method.
***********************************************************/
void solve_explicit_sequential(SimData_t *simData, 
                               int        xId);

/***********************************************************
* solve_implicit_sequential()
*-----------------------------------------------------------
* Solve the equation system 
*   A x = b
* using an implicit method.
***********************************************************/
void solve_implicit_sequential(SimData_t *simData, 
                               computeAx  cmpAx,
                               int        xId);

#endif /* SOLVER_LINEARSOLVER_H */
