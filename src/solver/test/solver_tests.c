#include <assert.h>

#include "aux/dbg.h"
#include "aux/minunit.h"
#include "solver/solver.h"
#include "solver/simData.h"
#include "solver/quadData.h"
#include "solver/typedefs.h"
#include "solver/dataIO.h"
#include "solver/gradients.h"

#include "solver_tests.h"

#include <stdio.h>
#include <stdlib.h>


#define CUR_TEST_NAME "solver_tests"


/************************************************************
* User-defined flow variable initialization
************************************************************/
void init_function(QuadData_t *quadData)
{
  int i;

  octDouble *xc = quadData->centroid;

  octDouble c[3] = {0.5, 0.5, 0.5};
  octDouble w    = 0.15;
  octDouble d[P4EST_DIM];

  octDouble r2 = 0.0;

  for (i = 0; i < P4EST_DIM; i++)
  {
    d[i] = xc[i] - c[i];
    r2  += d[i] * d[i];
  }

  octDouble arg = -0.5 * r2 / w / w;

  //quadData->vars[IS]  = sin(2.0*M_PI*r2/w)*exp(arg);
  quadData->vars[IS]  = exp(arg);
  quadData->vars[IVX] = 1.0;
  quadData->vars[IVY] = 1.0;
  quadData->vars[IVZ] = 0.0;
  quadData->vars[IRHO] = 1.0;

} /* init_function() */

/************************************************************
* Refinement function
************************************************************/
int refine_fn(p4est_t *p4est,
              p4est_topidx_t which_tree,
              p4est_quadrant_t *q)
{
  return 0;

  SimData_t     *simData  = (SimData_t*) p4est->user_pointer;
  SolverParam_t *solverParam = simData->solverParam;

  QuadData_t     *quadData = (QuadData_t *) q->p.user_data;
  octDouble      *xc       = quadData->centroid;

  if (q->level == solverParam->maxRefLvl)
    return 0;

  if (fabs(xc[0]-0.5) < 0.15 && fabs(xc[1]-0.5) < 0.15)
    return 1;

  return 0;
}

/************************************************************
* Coarsening function
************************************************************/
int coarse_fn(p4est_t *p4est,
              p4est_topidx_t which_tree,
              p4est_quadrant_t *children[])
{
  return 0;

  QuadData_t     *quadData;
  octDouble      *xc;
  p4est_quadrant_t *q;

  int i;

  for (i = 0; i < P4EST_CHILDREN; i++)
  {
    q        = children[i];
    quadData = (QuadData_t *) q->p.user_data;
    xc       = quadData->centroid;

    if (fabs(xc[0]-0.5) >= 0.25 && fabs(xc[1]-0.5) >= 0.25)
    {
      return 1;
    }
  }

  return 0;
}


/************************************************************
* Definition of unit test functions 
************************************************************/
char *test_solver_p4est(int argc, char *argv[])
{
  int                   refine_level = 0;
  int                   mpiret;
  int                   balance;
  int                   level;
  sc_MPI_Comm           mpicomm;
  p4est_t              *p4est;
  p4est_connectivity_t *conn;

  /*--------------------------------------------------------
  | Init MPI
  --------------------------------------------------------*/
  mpiret = sc_MPI_Init(&argc, &argv);
  SC_CHECK_MPI(mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /*--------------------------------------------------------
  | Optional functions for logging to processor zero
  --------------------------------------------------------*/
  sc_init(mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init(NULL, SC_LP_PRODUCTION);
  //P4EST_GLOBAL_PRODUCTIONF("UNIT TEST ",
  //   P4EST_DIM, P4EST_STRING);

  /*--------------------------------------------------------
  | Create a periodic connectivity
  --------------------------------------------------------*/
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_periodic ();
#else
  conn = p8est_connectivity_new_periodic ();
#endif

  /*--------------------------------------------------------
  | Create a forest that is not refined
  --------------------------------------------------------*/
  p4est = p4est_new(mpicomm, 
                    conn, 
                    0, 
                    NULL,
                    NULL);

  /*--------------------------------------------------------
  | Refine iteratively, load balancing at each iteration
  | --> important when starting wit unrefined forest
  --------------------------------------------------------*/
  for (level = 0; level < refine_level; ++level) {
    p4est_refine(p4est, 0, refine_fn, NULL);
    p4est_partition(p4est, 0, NULL);
  }

  /*--------------------------------------------------------
  | Call 2:1 balance to ensure that neighbors do not 
  | differ in size by more than a factor of 2.
  --------------------------------------------------------*/
  balance = 1;
  if (balance) {
    p4est_balance(p4est, P4EST_CONNECT_FACE, NULL);
    p4est_partition(p4est, 0, NULL);
  }

  /*--------------------------------------------------------
  | Destroy 
  --------------------------------------------------------*/
  p4est_destroy(p4est);
  p4est_connectivity_destroy(conn);

  /*--------------------------------------------------------
  | Verify memory leaks
  --------------------------------------------------------*/
  sc_finalize();

  /*--------------------------------------------------------
  | MPI ending
  --------------------------------------------------------*/
  mpiret = sc_MPI_Finalize();
  SC_CHECK_MPI(mpiret);

  return NULL;

} /* test_solver_p4est() */



/************************************************************
* Function to test the initialization and destruction of 
* the solver structures
************************************************************/
char *test_solver_init_destroy(int argc, char *argv[])
{

  SimData_t *simData = init_simData(argc, argv, 
                                    init_function,
                                    refine_fn,
                                    coarse_fn);

  writeSolutionVtk(simData, 0);

  solverRun(simData);

  destroy_simData(simData);

  return NULL;

} /* test_solver_init_destroy() */



