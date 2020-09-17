#-------------------------------------------------------------
#
# Findp4est
# ---------
#
# Options:
#
#   P4EST_ROOT           : Path to the P4EST base directory
#   P4EST_USE_STATIC     : Use static linking
#   P4EST_MPI_SUPPORT    : Enable MPI support
#
# This module defines the following variables:
#
#   P4EST_FOUND            : True if P4EST_INCLUDE_DIR is found
#   P4EST_INCLUDE_DIR      : where to find p4est.h, etc.
#   P4EST_INCLUDE_DIRS     : set when P4EST_INCLUDE_DIR found
#   P4EST_LIBRARIES        : the p4est libraries to link against.
#
# Authors:
#
#   - Tobias Dittmann <research@man-behind-moon.com>
#
#-------------------------------------------------------------
# save standard suffix
set(_P4EST_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

# share or static library
if( P4EST_USE_STATIC STREQUAL "ON" )
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  message(STATUS "Try to find static Sundials libraries.")
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
  message(STATUS "Try to find dynamic Sundials libraries.")
endif()

#-------------------------------------------------------------
# Find Root Path
#-------------------------------------------------------------
# find P4EST path
find_path(P4EST_ROOT_DIR
  include/p4est.h
  PATHS
    ${P4EST_ROOT}
#   $ENV{P4EST_HOME}
#   $ENV{P4ESTDIR}
)
message("-- p4est root directory: " ${P4EST_ROOT_DIR})

if( P4EST_ROOT_DIR )
  set( P4EST_ROOT_DIR ${P4EST_ROOT_DIR} CACHE STRING "The P4EST root directory." )
endif()


#-------------------------------------------------------------
# Find Include Paths
#-------------------------------------------------------------
find_path(P4EST_INCLUDE_DIR
  p4est.h
  PATHS
    ${P4EST_ROOT_DIR}/include
#   ${INCLUDE_INSTALL_DIR}
)
message("-- p4est include directory: " ${P4EST_INCLUDE_DIR})


#-------------------------------------------------------------
# Find Libraries
#-------------------------------------------------------------
# SC auxiliary library
find_library(P4EST_SC_LIBRARY
  sc
  PATHS
    ${P4EST_ROOT_DIR}/lib/
)
message("-- p4est SC library: " ${P4EST_SC_LIBRARY})

# p4est library
find_library(P4EST_LIBRARY
  p4est
  PATHS
    ${P4EST_ROOT_DIR}/lib/
)
message("-- p4est library:    " ${P4EST_LIBRARY})

# copy libraries
set( P4EST_LIBRARIES "${P4EST_SC_LIBRARY} ${P4EST_LIBRARY}" )
message("-- p4est libraries:  " ${P4EST_LIBRARIES})

mark_as_advanced(
  P4EST_ROOT_DIR
  P4EST_INCLUDE_DIR
  P4EST_SC_LIBRARY
  P4EST_LIBRARY
)

# restore standard suffix
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_P4EST_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

