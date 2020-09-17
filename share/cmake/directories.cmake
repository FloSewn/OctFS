set( SRC         ${CMAKE_SOURCE_DIR}/src         )
set( INC         ${CMAKE_SOURCE_DIR}/include     )
set( LIB         ${CMAKE_SOURCE_DIR}/lib         )
set( BIN         ${CMAKE_SOURCE_DIR}/bin         )
set( CMAKE_SHARE ${CMAKE_SOURCE_DIR}/share/cmake )

# Define relative module path
string( LENGTH ${CMAKE_SOURCE_DIR} CMAKE_SOURCE_DIR_LEN )

# dependecies
set( P4EST_HOME /home/florian/Programs/p4est-2.2/local )
set( MPI_HOME   /usr/lib/x86_64-linux-gnu/openmpi      )
