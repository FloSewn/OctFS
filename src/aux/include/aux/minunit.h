#undef NDEBUG
#ifndef _minunit_h
#define _minunit_h

#include <stdio.h>
#include <stdlib.h>

#include "aux/dbg.h"

int tests_run;

#define mu_suite_start() char *message = NULL

#define mu_assert(test, message) if (!(test)) {\
  log_err(message); return message; }

#define mu_run_test(test, argc, argv) debug("\n-----%s", " " #test); \
  message = test(argc, argv); tests_run++; if (message) return message;

#define mu_print_tests_run() debug("\nTESTS RUN: %d\n", tests_run)

#define RUN_TESTS(name) int main(int argc, char *argv[]) {\
  argc = 0;argc++;\
  debug("----- RUNNING: %s", argv[0]);\
  printf("----\nRUNNING: %s\n", argv[0]);\
  char *result = name();\
  if (result != 0) {\
    printf("FAILED: %s\n", result);\
  }\
  else {\
    printf("ALL TESTS PASSED\n");\
  }\
  debug("Tests run: %d\n", tests_run);\
  exit(result != 0);\
}


#endif

//argc = 1; 
