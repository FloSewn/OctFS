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
#ifndef SOLVER_UTIL_H
#define SOLVER_UTIL_H

#include <math.h>
#include <stdio.h>

#include "solver/typedefs.h"

/***********************************************************
* True / False conditionals
***********************************************************/
#ifndef TRUE 
#define TRUE 1
#endif

#ifndef FALSE 
#define FALSE 0
#endif

/***********************************************************
* Own print function
***********************************************************/
#define octPrint(M, ...) fprintf(stdout, "> " M "\n",\
    ##__VA_ARGS__)

/***********************************************************
* Define maximum and minimum floats / doubles to 
* prevent overflows
***********************************************************/
#ifndef FLT_MIN
#define FLT_MIN 1.e-20
#endif

#ifndef FLT_MAX
#define FLT_MAX 1.e+20
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.e-7
#endif

#ifndef DBL_MIN
#define DBL_MIN 1.e-200
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.e+200
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.e-16
#endif

/***********************************************************
* Set final constants depending on precision
* -> Defined by OCT_USE_DOUBLE flag in oct_typedefs.h
***********************************************************/
#ifndef SYSTEM_MAXNUM
#ifdef OCT_USE_DOUBLE
#define SYSTEM_MAXNUM DBL_MAX
#else
#define SYSTEM_MAXNUM FLT_MAX
#endif
#endif

#ifndef SYSTEM_MINNUM
#ifdef OCT_USE_DOUBLE
#define SYSTEM_MINNUM DBL_MIN
#else
#define SYSTEM_MINNUM FLT_MIN
#endif
#endif

#ifndef SYSTEM_EPSILON
#ifdef OCT_USE_DOUBLE
#define SYSTEM_EPSILON DBL_EPSILON
#else
#define SYSTEM_EPSILON FLT_EPSILON
#endif
#endif

/***********************************************************
*
***********************************************************/
#ifndef SQR
#define SQR(x) ( (x) * (x) ) 
#endif

#ifndef POW3
#define POW3(x) ( (x) * (x) * (x) ) 
#endif


#ifndef MAX
#define MAX(a, b) ( (a) > (b) ? (a) : (b) ) 
#endif

#ifndef MIN
#define MIN(a, b) ( (a) < (b) ? (a) : (b) ) 
#endif


#ifndef MAX0
#define MAX0(a) ( (a) > 0 ? (a) : 0 ) 
#endif

#ifndef MIN0
#define MIN0(a) ( (a) < 0 ? (a) : 0 ) 
#endif


#ifndef ABS
#define ABS(a) ( (a) < 0 ? -1 * (a) : (a) ) 
#endif



#endif /* SOLVER_UTIL_H */
