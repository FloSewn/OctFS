#-------------------------------------------------------------
#
# Findmpi
# ---------
#
# Options:
#
#   MPI_ROOT           : Path to the P4EST base directory
#   MPI_USE_STATIC     : Use static linking
#
# This module defines the following variables:
#
#   MPI_FOUND            : True if MPI_INCLUDE_DIR is found
#   MPI_INCLUDE_DIR      : where to find mpi.h, etc.
#   MPI_INCLUDE_DIRS     : set when MPI_INCLUDE_DIR found
#   MPI_LIBRARIES        : the mpi libraries to link against.
#
# Authors:
#
#   - Florian Setzwein 
#
#-------------------------------------------------------------
# save standard suffix
set(_MPI_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

# share or static library
if( MPI_USE_STATIC STREQUAL "ON" )
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  message(STATUS "Try to find static Sundials libraries.")
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
  message(STATUS "Try to find dynamic Sundials libraries.")
endif()

#-------------------------------------------------------------
# Find Root Path
#-------------------------------------------------------------
# find mpi path
find_path(MPI_ROOT_DIR
  include/mpi.h
  PATHS
    ${MPI_ROOT}
)
message("-- mpi root directory: " ${MPI_ROOT_DIR})

#-------------------------------------------------------------
# Find Include Paths
#-------------------------------------------------------------
find_path(MPI_INCLUDE_DIR
  mpi.h
  PATHS
    ${MPI_ROOT_DIR}/include
)
message("-- mpi include directory: " ${MPI_INCLUDE_DIR})

#-------------------------------------------------------------
# Find Libraries
#-------------------------------------------------------------
find_library(MPI_LIBRARY
  mpi
  PATHS
    ${MPI_ROOT_DIR}/lib/
)
message("-- mpi library:    " ${MPI_LIBRARY})

# copy libraries
set( MPI_LIBRARIES "${MPI_LIBRARY}" )
message("-- mpi libraries:  " ${MPI_LIBRARIES})

mark_as_advanced(
  MPI_ROOT_DIR
  MPI_INCLUDE_DIR
  MPI_LIBRARY
)

# restore standard suffix
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MPI_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

