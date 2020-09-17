#############################################################
# Define dependencies                                       #
#############################################################

# extend CMake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SHARE}/modules/")

# Mandatory libraries
set( LIB_P4EST_ROOT ${P4EST_HOME} CACHE PATH "Defines the p4est root directory." )
set( P4EST_ROOT ${P4EST_HOME} )
message("-- p4est root:  " ${P4EST_ROOT} )

set( LIB_MPI_ROOT ${MPI_HOME} CACHE PATH "Defines the mpi root directory." )
set( MPI_ROOT ${MPI_HOME} )
message("-- mpi root:  " ${MPI_ROOT} )

# configure module
set( P4EST_USE_STATIC ON )
set( MPI_USE_STATIC OFF )

# find module
find_package( p4est REQUIRED )
find_package( mpi REQUIRED )
