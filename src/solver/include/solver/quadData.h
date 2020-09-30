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
#ifndef SOLVER_QUADDATA_H
#define SOLVER_QUADDATA_H

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
* Structure containing properties for every quad
*   > Accessed through q->p.user_data
***********************************************************/
typedef struct QuadData_t
{
  /*--------------------------------------------------------
  | Quad geometry data
  --------------------------------------------------------*/
  // Vertices 
#ifdef P4_TO_P8
  octDouble xyz[8][3];
#else
  octDouble xyz[4][2];
#endif
  // Octahedron centroid 
  octDouble centroid[P4EST_DIM];
  // Octahedron volume 
  octDouble volume;
  // Octahedron face normals 
  octDouble normals[2*P4EST_DIM][P4EST_DIM];
  // Octahedron face centroids 
  octDouble face_centroids[2*P4EST_DIM][P4EST_DIM];

  /*--------------------------------------------------------
  | Quad flow data
  --------------------------------------------------------*/
  // Masfluxes
  octDouble mflux[2*P4EST_DIM];
  // State variables 
  octDouble vars[OCT_MAX_VARS];
  // State variable gradients 
  octDouble grad_vars[OCT_MAX_VARS][P4EST_DIM];
  // State variable  buffers 
  octDouble vars_buf[OCT_MAX_VARS];

  /*--------------------------------------------------------
  | Solver buffer data
  --------------------------------------------------------*/
  octDouble *Ax_p;
  octDouble  Ax[OCT_MAX_VARS];
  octDouble  b[OCT_MAX_VARS];
  octDouble  res[OCT_MAX_VARS];

} QuadData_t;

/***********************************************************
* init_quadData()
*-----------------------------------------------------------
* Initializes the quad data structure
***********************************************************/
void init_quadData(p4est_t *p4est,
                   p4est_topidx_t which_tree,
                   p4est_quadrant_t *q);

/***********************************************************
* init_quadFlowData()
*-----------------------------------------------------------
* Initializes the quad flow data structure
***********************************************************/
void init_quadFlowData(QuadData_t *quadData);

/***********************************************************
* init_quadGeomData3d()
*-----------------------------------------------------------
* Initializes the 3D quad geometry data structure
***********************************************************/
void init_quadGeomData3d(p4est_t          *p4est,
                         p4est_topidx_t    which_tree,
                         p4est_quadrant_t *q, 
                         QuadData_t       *quadData);

/***********************************************************
* init_quadGeomData2d()
*-----------------------------------------------------------
* Initializes the 2D quad geometry data structure
***********************************************************/
void init_quadGeomData2d(p4est_t          *p4est,
                         p4est_topidx_t    which_tree,
                         p4est_quadrant_t *q, 
                         QuadData_t       *quadData);

/***********************************************************
* interpQuadData()
*-----------------------------------------------------------
* Function for initialzing the state variables of incoming
* quadrants from outgoing quadrants. 
*
* The functions p4est_refine_ext(), p4est_coarsen_ext() and
* p4est_balance_ext() take as an arguments a p4est_replace_t
* callback function.
* This function allows to setup the quadrant data of
* incoming quadrants from the data of outgoing quadrants, 
* before the outgoing data is destroyed.
* This function matches the p4est_replace_t prototype.
*
* In this example, we linearly interpolate the state 
* variable of a quadrant that is refined to its children
* and we average the children midpoints that are being 
* coarsened to the paren
*-----------------------------------------------------------
* Arguments:
* *p4est        : the forest
*  which_tree   : the tree in the forest containing a 
*                 child
*  num_outgoing : the number of quadrants that are being
*                 replaced:
*                 Either 1 if a quadrant is being refined 
*                 or P4EST_CHILDREN if a family of 
*                 children are being coarsened
*  outgoing     : the outgoing quadrants
*  num_incoming : the number of quadrants that are being
*                 added:
*                 Either P4EST_CHILDREN if a quadrant is 
*                 being refined, or 1 if a family of 
*                 children are being coarsened.
*  incoming     : quadrants whose data is initialized.
*
***********************************************************/
void interpQuadData(p4est_t          *p4est, 
                    p4est_topidx_t    which_tree,
                    int               num_outgoing,
                    p4est_quadrant_t *outgoing[],
                    int               num_incoming, 
                    p4est_quadrant_t *incoming[]);


#endif /* SOLVER_QUADDATA_H */
