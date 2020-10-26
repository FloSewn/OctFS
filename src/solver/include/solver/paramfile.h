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
#ifndef AUX_PARAMFILE_H
#define AUX_PARAMFILE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "aux/bstrlib.h"
#include "solver/typedefs.h"

#define FILIO_ERR -1
#define OCT_MAX_PARAMETERS 100
#define OCT_MAX_PARAM_LENGTH 50


struct octParam;
struct octParamInst;

/***********************************************************
* Parameter type identifier
***********************************************************/
typedef enum
{
  INTVAL, /* single integer value */
  DBLVAL, /* single double value */
  STRVAL, /* single string value */
  DBLVEC  /* vector of doubles */
} paramType;

/*************************************************************
* Parameter file instruction
*************************************************************/
typedef struct octParamInst {

  /* Instruction string  */
  char       inst[OCT_MAX_PARAM_LENGTH];      

  /* parameter value */
  void      *value;

  /* Parameter type */
  paramType  pType;     

  /* Flag if parameter is mandatory */
  octBool    mandatory;

  /* default values */
  int        intDefault;
  octDouble  dblDefault;
  char      *strDefault;

} octParamInst;

/*************************************************************
* Parameter file structure
*************************************************************/
typedef struct octParam {
  const char      *path;    /* Path of file                 */
  bstring          txt;     /* bstring with file data       */
  struct bstrList *txtlist; /* file, splitted for newlines  */

  long             length;  /* Number of chars in total file*/
                            /* -> including '\0' at end     */
  int              nlines;  /* Number of lines in total file*/

} octParam;


/*************************************************************
* octParam_readParamfile()
*-------------------------------------------------------------
* Function to read the parameter file and initialize  
* respective parameters
*************************************************************/
int octParam_readParamfile(SimData_t  *simData,
                           const char *filePath);

/*************************************************************
* octParam_initParameters()
*-------------------------------------------------------------
* Function to read the parameter file and initialize  
* respective parameters
*************************************************************/
int octParam_initParameters(SimData_t *simData,
                            octParam  *paramFile);

/*************************************************************
* Function to create a new parameter file reader structure
*************************************************************/
octParam *octParam_create(const char *file_path);

/*************************************************************
* Function to destroy a file reader structure
*************************************************************/
int octParam_destroy(octParam *file);

/*************************************************************
* Function returns a bstring list of lines, that 
* do not contain a certain specifier
*************************************************************/
struct bstrList *octParam_popLinesWith(struct bstrList *txtlist,
                                      const char *fltr);

/*************************************************************
* Function returns a bstring list of lines, that 
* do contain a certain specifier
*************************************************************/
struct bstrList *octParam_getLinesWith(struct bstrList *txtlist,
                                      const char *fltr);

/*************************************************************
* Function searches for a specifier <fltr> in a bstrList.
* The parameter behind the specifier is then extracted 
* from the file and stored into <value>.
* The value is casted to a prescribed type
* type = 0: integer
* type = 1: double
* type = 2: string
*
* Returns 0 if specifier was not found in the file.
* Otherwise, it returns the number of times, the 
* specifier was found.
* Returns -1 on errors.
*************************************************************/
int octParam_extractParam(struct bstrList *txtlist,
                         const char *fltr, int type,
                         void *value);

/*************************************************************
* Function searches for a specifier <fltr> in a bstrList.
* The string behind the specifier is then extracted 
* from the file and processed as an array of values 
* and stored in <value>.
* The values are casted to a prescribed type
* type = 0: integer
* type = 1: double
* type = 2: string
*
* Returns 0 if specifier was not found in the file.
* Otherwise, it returns the number of times, the 
* specifier was found.
* Returns -1 on errors.
*************************************************************/
int octParam_extractArray(struct bstrList *txtlist,
                         const char *fltr, int type,
                         void *value);

#endif
