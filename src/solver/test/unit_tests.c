#include <assert.h>

#include "aux/dbg.h"
#include "aux/minunit.h"

#include "solver_tests.h"


/************************************************************
* Run all unit test functions
************************************************************/
char *all_tests(int argc, char *argv[])
{
  mu_suite_start();

  mu_run_test(test_solver_init_destroy, argc, argv);

  return NULL;
}

/************************************************************
* Main function to run unit tests
************************************************************/
int main(int argc, char *argv[])
{
  debug("----- RUNNING %s\n", argv[0]);

  char *result;
  result = all_tests(argc, argv);

  if (result != NULL)
  {
    debug("\nTESTS FAILED!\n");
  }
  else
  {
    debug("\nALL TESTS PASSED!\n");
  }

  mu_print_tests_run();

  return 0;

}
