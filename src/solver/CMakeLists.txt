##############################################################
# MODULE: SOLVER
##############################################################
# Define directories
set( SOLVER_DIR    ${SRC}/solver         )
set( SOLVER_INC    ${SOLVER_DIR}/include )
set( SOLVER_SRC    ${SOLVER_DIR}/src     )

# Define source files
set( SOLVER_MAIN
  ${SOLVER_SRC}/simData.c
  ${SOLVER_SRC}/quadData.c
  ${SOLVER_SRC}/solver.c
  ${SOLVER_SRC}/refine.c
  ${SOLVER_SRC}/coarsen.c
  ${SOLVER_SRC}/dataIO.c
  ${SOLVER_SRC}/gradients.c
  )

##############################################################
# LIBRARY: SOLVER
##############################################################
# define name
set( SOLVER_LIB solver )

# add sources to library
add_library( ${SOLVER_LIB} STATIC ${SOLVER_MAIN} )

# set includes
target_include_directories( ${SOLVER_LIB}
  PUBLIC $<BUILD_INTERFACE:${SOLVER_INC}>
  PRIVATE ${SOLVER_SRC}
  PUBLIC ${P4EST_INCLUDE_DIR}
  PUBLIC ${MPI_INCLUDE_DIR}
)

# set libraries
target_link_libraries( ${SOLVER_LIB}
  PUBLIC aux
  INTERFACE m
  INTERFACE lapack
  INTERFACE cblas
  INTERFACE f77blas
  INTERFACE atlas
  INTERFACE z
  PUBLIC ${P4EST_LIBRARY}
  PUBLIC ${P4EST_SC_LIBRARY}
  PUBLIC ${MPI_LIBRARY}
)

install( TARGETS solver DESTINATION ${LIB} )


##############################################################
# TESTS: solver
##############################################################
set( TESTEXE_SOLVER solver_test        )
set( TESTDIR_SOLVER ${SOLVER_DIR}/test )

add_executable( ${TESTEXE_SOLVER}
  ${TESTDIR_SOLVER}/solver_tests.c
  ${TESTDIR_SOLVER}/unit_tests.c
)

target_include_directories( ${TESTEXE_SOLVER}
  PUBLIC $<BUILD_INTERFACE:${INC_SOLVER}>
  PUBLIC ${P4EST_INCLUDE_DIR}
  PUBLIC ${MPI_INCLUDE_DIR}
)

target_link_libraries( ${TESTEXE_SOLVER}
  PUBLIC solver
  PUBLIC aux
  INTERFACE m
  INTERFACE lapack
  INTERFACE cblas
  INTERFACE f77blas
  INTERFACE atlas
  INTERFACE z
  PUBLIC ${P4EST_LIBRARY}
  PUBLIC ${P4EST_SC_LIBRARY}
  PUBLIC ${MPI_LIBRARY}
)


# Install executables
install( TARGETS ${TESTEXE_SOLVER} RUNTIME DESTINATION ${BIN} )



