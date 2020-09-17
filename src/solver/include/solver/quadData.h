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
* Structure containing geometric properties for every quad
***********************************************************/
typedef struct QuadGeomData_t 
{
  /* Vertices */
#ifdef P4_TO_P8
  octDouble xyz[8][3];
#else
  octDouble xyz[4][2];
#endif

  /* Octahedron centroid */
  octDouble centroid[P4EST_DIM];
  /* Octahedron volume */
  octDouble volume;
  /* Octahedron face normals */
  octDouble normals[2*P4EST_DIM][P4EST_DIM];
  /* Octahedron face centroids */
  octDouble face_centroids[2*P4EST_DIM][P4EST_DIM];

} QuadGeomData_t;


/***********************************************************
* Structure containing the flow field variable data 
* for every quad
***********************************************************/
typedef struct QuadFlowData_t
{
  /* State variables */
  octDouble vars[OCT_MAX_VARS];

  /* State variable gradients */
  octDouble grad_vars[OCT_MAX_VARS][P4EST_DIM];

  /* State variable  buffers */
  octDouble vars_buf[OCT_MAX_VARS];

} QuadFlowData_t;


/***********************************************************
* Structure containing properties for every quad
*   > Accessed through q->p.user_data
***********************************************************/
typedef struct QuadData_t
{
  QuadGeomData_t geomData;
  QuadFlowData_t flowData;

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
void init_quadFlowData(p4est_t *p4est,
                       p4est_topidx_t which_tree,
                       p4est_quadrant_t *q, 
                       QuadFlowData_t *flowData);

/***********************************************************
* init_quadGeomData3d()
*-----------------------------------------------------------
* Initializes the 3D quad geometry data structure
***********************************************************/
void init_quadGeomData3d(p4est_t          *p4est,
                         p4est_topidx_t    which_tree,
                         p4est_quadrant_t *q, 
                         QuadGeomData_t   *geomData);

/***********************************************************
* init_quadGeomData2d()
*-----------------------------------------------------------
* Initializes the 2D quad geometry data structure
***********************************************************/
void init_quadGeomData2d(p4est_t          *p4est,
                         p4est_topidx_t    which_tree,
                         p4est_quadrant_t *q, 
                         QuadGeomData_t   *geomData);


#endif /* SOLVER_QUADDATA_H */
