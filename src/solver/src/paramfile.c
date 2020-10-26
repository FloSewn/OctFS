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
#include "aux/dbg.h"
#include "aux/bstrlib.h"
#include "solver/paramfile.h"
#include "solver/simData.h"

/*************************************************************
* octParam_readParamfile()
*-------------------------------------------------------------
* Function to read the parameter file and initialize  
* respective parameters
*************************************************************/
int octParam_readParamfile(SimData_t  *simData,
                           const char *filePath)
{
  /*----------------------------------------------------------
  | Load parameter file and clear all comments
  ----------------------------------------------------------*/
  octParam *paramFile = octParam_create(filePath);
  struct bstrList *buf = octParam_popLinesWith(paramFile->txtlist, "#");
  bstrListDestroy(paramFile->txtlist);
  paramFile->txtlist = buf;

  /*----------------------------------------------------------
  |
  ----------------------------------------------------------*/
  int stopSim = FALSE;
  stopSim = octParam_initParameters(simData, paramFile);

  /*----------------------------------------------------------
  | Free memory
  ----------------------------------------------------------*/
  octParam_destroy( paramFile );

  return stopSim;

} /* octParam_readParamfile() */

/*************************************************************
* octParam_initParameters()
*-------------------------------------------------------------
* Function to read the parameter file and initialize  
* respective parameters
*************************************************************/
int octParam_initParameters(SimData_t *simData,
                            octParam  *paramFile)
{
  SimParam_t    *simParam    = simData->simParam;
  SolverParam_t *solverParam = simData->solverParam;

  /*----------------------------------------------------------
  | Define simulation parameter instructions
  ----------------------------------------------------------*/
  octParamInst simParamInst[OCT_MAX_PARAMETERS] = 
  {
    {"Simulation time step [s]:", 
     &simParam->timestep, DBLVAL, TRUE, 
     -1, -1.0, NULL},
    {"Total simulation time [s]:",
     &simParam->simTimeTot, DBLVAL, TRUE, 
     -1, -1.0, NULL},
    {"Temporal discretization scheme:",
     &simParam->tempScheme, STRVAL, FALSE, 
     -1, -1.0, "Euler-Forward"},
    {"Reference kinematic viscosity [Pa*s]:",
     &simParam->viscosity, DBLVAL, FALSE, 
     -1, 1.0E-5, NULL},
  };

  /*----------------------------------------------------------
  | Define solver parameter instructions
  ----------------------------------------------------------*/


  
  /*----------------------------------------------------------
  | Read simulation parameter instructions
  ----------------------------------------------------------*/
  int i, nvals;
  octBool stopSim = FALSE;

  for (i = 0; i < OCT_MAX_PARAMETERS; i++)
  {
    if (strlen(simParamInst[i].inst) == 0)
      continue;

    nvals = octParam_extractParam(paramFile->txtlist,
                                  simParamInst[i].inst,
                                  simParamInst[i].pType,
                                  simParamInst[i].value);

    /*--------------------------------------------------------
    | Handle missing parameters
    --------------------------------------------------------*/
    if (nvals < 1 && simParamInst[i].mandatory == TRUE)
    {
      octPrint("[ERROR]: MISSING PARAMETER");
      octPrint("%s <UNDEFINED>", simParamInst[i].inst);
      stopSim = TRUE;
    }
    else if (nvals < 1 && simParamInst[i].mandatory == FALSE)
    {
      if (simParamInst[i].pType == INTVAL)
      {
        *(int*)simParamInst[i].value = simParamInst[i].intDefault;
      }
      else if (simParamInst[i].pType == DBLVAL)
      {
        *(octDouble*)simParamInst[i].value = simParamInst[i].dblDefault;
      }
      else if (simParamInst[i].pType == STRVAL)
      {
        *(char*)simParamInst[i].value = simParamInst[i].strDefault;
      }
    }

    /*--------------------------------------------------------
    | Output parameters to user
    --------------------------------------------------------*/
    if (simParamInst[i].pType == INTVAL)
    {
      int *printVal = (int*)simParamInst[i].value;
      octPrint("%s %d", simParamInst[i].inst, *printVal);
    }
    else if (simParamInst[i].pType == DBLVAL)
    {
      octDouble *printVal = (octDouble*)simParamInst[i].value;
      octPrint("%s %e", simParamInst[i].inst, *printVal);
    }
    else if (simParamInst[i].pType == STRVAL)
    {
      char *printVal = (char*)(*(bstring*)simParamInst[i].value)->data;
      octPrint("%s %s", simParamInst[i].inst, printVal);
    }

  }


  return stopSim;

} /* octParam_initParameters() */

/*************************************************************
* Function to create a new parameter file reader
*
* Reference:
* https://stackoverflow.com/questions/14002954/c-programming
* -how-to-read-the-whole-file-contents-into-a-buffer
*************************************************************/
octParam *octParam_create(const char *file_path)
{
  /*---------------------------------------------------------
  | Allocate memory for txtio structure 
  ---------------------------------------------------------*/
  octParam *txtfile = calloc(1, sizeof(octParam));
  check_mem(txtfile);

  txtfile->path = file_path;

  /*---------------------------------------------------------
  | Open text file and copy its data 
  ---------------------------------------------------------*/
  FILE *fptr = NULL;
  fptr = fopen(txtfile->path, "rb");
  check(fptr, "Failed to open %s.", txtfile->path);

  /* Estimate length of chars in whole file                */
  fseek(fptr, 0, SEEK_END);
  long length = ftell(fptr);
  fseek(fptr, 0, SEEK_SET);

  /* Read total file into buffer                           */
  char *buffer = (char *) malloc(length + 1);
  buffer[length] = '\0';
  fread(buffer, 1, length, fptr);

  /* Copy relevant data to octParam structure               */
  bstring bbuffer = bfromcstr( buffer );
  txtfile->txt    = bbuffer;
  txtfile->length = length + 1;

  /* Split buffer according to newlines                    */
  char splitter = '\n';
  txtfile->txtlist = bsplit(bbuffer, splitter);
  txtfile->nlines = txtfile->txtlist->qty;

  fclose(fptr);
  free(buffer);

  return txtfile;
error:
  return NULL;
}

/*************************************************************
* Function to destroy a file reader structure
*************************************************************/
int octParam_destroy(octParam *file)
{
  bstrListDestroy(file->txtlist);
  bdestroy(file->txt);
  free(file);
  return 0;
}

/*************************************************************
* Function returns a bstring list of lines, that 
* do not contain a specifier
*************************************************************/
struct bstrList *octParam_popLinesWith(struct bstrList *txtlist, 
                                      const char *fltr)
{
  bstring *fl_ptr = txtlist->entry;
  bstring  bfltr  = bfromcstr(fltr); 

  bstring *bl_ptr     = NULL;
  struct bstrList *bl = NULL;
  int i;
  int hits = 0;
  int cnt  = 0;
  int *marker = (int*) calloc(txtlist->qty, sizeof(int));
  check_mem(marker);

  /*----------------------------------------------------------
  | Fill array of markers with line numbers, that do not 
  | contain the filter string
  ----------------------------------------------------------*/
  for (i = 0; i < txtlist->qty; i++) 
  {
    if ( binstr(fl_ptr[i], 0, bfltr) == BSTR_ERR )
    {
      hits += 1;
      marker[i] = i;
    }
    else
    {
      marker[i] = -1;
    }
  }

  /*----------------------------------------------------------
  | Create new bstrList that will contained all filtered lines
  ----------------------------------------------------------*/
  bl = bstrListCreate();
  bstrListAlloc(bl, hits);

  /*----------------------------------------------------------
  | Copy marked lines into new bstrList
  ----------------------------------------------------------*/
  bl_ptr = bl->entry;

  for (i = 0; i < txtlist->qty; i++) 
  {
    if (marker[i] >= 0)
    {
      const int curline = marker[i];
      bl_ptr[cnt] = bstrcpy(fl_ptr[curline]);
      bl->qty += 1;
      cnt += 1;
    }
  }

  /*----------------------------------------------------------
  | Cleanup
  ----------------------------------------------------------*/
  bdestroy(bfltr);
  free(marker);

  return bl;
error:

  bdestroy(bfltr);
  free(marker);
  return NULL;
}

/*************************************************************
* Function returns a bstring list of lines, that 
* do contain a certain specifier
*************************************************************/
struct bstrList *octParam_getLinesWith(struct bstrList *txtlist,
                                      const char *fltr)
{
  bstring *fl_ptr = txtlist->entry;
  bstring  bfltr  = bfromcstr( fltr ); 

  bstring *bl_ptr     = NULL;
  struct bstrList *bl = NULL;

  int i;
  int hits = 0;
  int cnt  = 0;
  int *marker = (int*) calloc(txtlist->qty, sizeof(int));
  check_mem(marker);

  /*----------------------------------------------------------
  | Fill array of markers with line numbers, that  
  | contain the filter string
  ----------------------------------------------------------*/
  for (i = 0; i < txtlist->qty; i++) 
  {
    if ( binstr(fl_ptr[i], 0, bfltr) != BSTR_ERR )
    {
      hits += 1;
      marker[i] = i;
    }
    else
    {
      marker[i] = -1;
    }
  }

  /*----------------------------------------------------------
  | Create new bstrList that will contained all filtered lines
  ----------------------------------------------------------*/
  bl = bstrListCreate();
  bstrListAlloc(bl, hits);

  /*----------------------------------------------------------
  | Copy marked lines into new bstrList
  ----------------------------------------------------------*/
  bl_ptr = bl->entry;

  for (i = 0; i < txtlist->qty; i++) 
  {
    if (marker[i] >= 0)
    {
      const int curline = marker[i];
      bl_ptr[cnt] = bstrcpy(fl_ptr[curline]);
      bl->qty += 1;
      cnt += 1;
    }
  }

  /*----------------------------------------------------------
  | Cleanup
  ----------------------------------------------------------*/
  bdestroy(bfltr);
  free(marker);

  return bl;

error:
  bdestroy(bfltr);
  free(marker);
  return NULL;
}

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
                         void *value)
{
  int nfound = 0;
  bstring line, valstr;
  struct bstrList *fltTxt = NULL;
  struct bstrList *blstextr = NULL;
  bstring bfltr, bextr;

  /*----------------------------------------------------------
  | Get all lines, containing the specifier
  ----------------------------------------------------------*/
  fltTxt = octParam_getLinesWith(txtlist, fltr);
  nfound = fltTxt->qty;

  /*----------------------------------------------------------
  | Return if specifier is not found
  ----------------------------------------------------------*/
  if (fltTxt->qty < 1)
  {
    bstrListDestroy(fltTxt);
    return 0;
  }

  /*----------------------------------------------------------
  | Take last string, in which specifier was found
  ----------------------------------------------------------*/
  line  = fltTxt->entry[fltTxt->qty - 1];
  bfltr = bfromcstr( fltr ); 
  
  int off = binstr(line, 0, bfltr); 
  int len = bfltr->slen;

  bextr = bmidstr( line, off+len, line->slen );

  /*----------------------------------------------------------
  | Remove leading whitespaces and copy first value
  ----------------------------------------------------------*/
  valstr = bextr;

  if (type == INTVAL)
  {
    *(int*)value = atoi((char*) valstr->data);
  }
  else if (type == DBLVAL)
  {
    *(double*)value = atof((char*) valstr->data);
  }
  else if (type == STRVAL)
  {
    *(bstring*)value = bfromcstr((char*) valstr->data );
  }
  else
    log_err("Wrong type definition.");

  /*----------------------------------------------------------
  | Cleanup
  ----------------------------------------------------------*/
  bdestroy( bextr );
  bdestroy( bfltr );
  bstrListDestroy( fltTxt );
  bstrListDestroy( blstextr );

  return nfound;

}

/*************************************************************
* Function searches for a specifier <fltr> in a bstrList.
* The string behind the specifier is then extracted 
* from the file and processed as an array of values 
* and stored in <value>.
* The values are casted to a prescribed type
* type = 0: integer
* type = 1: double
* type = 2: bstrList (not working )
*
* Returns 0 if specifier was not found in the file.
* Otherwise, it returns the number of times, the 
* specifier was found.
* Returns -1 on errors.
*************************************************************/
int octParam_extractArray(struct bstrList *txtlist,
                         const char *fltr, int type,
                         void *value)
{

  int i;
  int nfound = 0;
  bstring line;
  struct bstrList *fltTxt = NULL;
  bstring bfltr, bextr;

  /*----------------------------------------------------------
  | Get all lines, containing the specifier
  ----------------------------------------------------------*/
  fltTxt = octParam_getLinesWith(txtlist, fltr);
  nfound = fltTxt->qty;

  /*----------------------------------------------------------
  | Return if specifier is not found
  ----------------------------------------------------------*/
  if (nfound < 1)
  {
    bstrListDestroy(fltTxt);
    return 0;
  }

  /*----------------------------------------------------------
  | Take last string, in which specifier was found
  ----------------------------------------------------------*/
  line  = fltTxt->entry[fltTxt->qty - 1];
  bfltr = bfromcstr( fltr ); 

  int off = binstr(line, 0, bfltr); 
  int len = bfltr->slen;

  bextr = bmidstr( line, off+len, line->slen );

  /*----------------------------------------------------------
  | Remove leading whitespaces and copy first value
  ----------------------------------------------------------*/
  bstring wsfnd = bfromcstr( " " );
  bstring wsrpl = bfromcstr( "" );
  bfindreplace(bextr, wsfnd, wsrpl, 0);

  /*----------------------------------------------------------
  | Split into list of string -> comma is separator
  ----------------------------------------------------------*/
  struct bstrList *arrStr = bsplit(bextr, ',');
  int nEntries            = arrStr->qty;
  bstring *arr_ptr        = arrStr->entry;

  if (type == 0)
  {
    int *array = calloc(nEntries, sizeof(int));

    for (i = 0; i < nEntries; i++) 
      array[i] = atoi((char*) arr_ptr[i]->data);

    *(int**)value = array;

  }
  else if (type == 1)
  {
    double *array = calloc(nEntries, sizeof(double));

    for (i = 0; i < nEntries; i++) 
      array[i] = atof((char*) arr_ptr[i]->data);

    *(double**)value = array;

  }
  else if (type == 2)
  {
    *(struct bstrList**)value = arrStr;
  }
  else
    log_err("Wrong type definition.");

  /*----------------------------------------------------------
  | Cleanup
  ----------------------------------------------------------*/
  bdestroy( wsfnd );
  bdestroy( wsrpl );
  bdestroy( bextr );
  bdestroy( bfltr );
  bstrListDestroy( fltTxt );
  bstrListDestroy( arrStr);

  return nfound;


} /* octParam_extractArray() */

