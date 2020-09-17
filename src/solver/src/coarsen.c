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
#include "solver/typedefs.h"
#include "solver/util.h"
#include "solver/coarsen.h"
#include "solver/quadData.h"
#include "solver/simData.h"


/***********************************************************
* globalCoarsening()
*-----------------------------------------------------------
* This is the general function to control the coarsening
* of trees.
***********************************************************/
int globalCoarsening(p4est_t          *p4est,
                     p4est_topidx_t    which_tree,
                     p4est_quadrant_t *children[])
{
  SimData_t  *simData  = (SimData_t*) p4est->user_pointer;
  SimParam_t *simParam = simData->simParam;

  int coarsen = 0;

  if (simParam->usrCoarseFun != NULL)
  {
    coarsen |= simParam->usrCoarseFun(p4est, 
                                      which_tree, 
                                      children);
  }



  return coarsen;

} /* globalCoarsening() */
