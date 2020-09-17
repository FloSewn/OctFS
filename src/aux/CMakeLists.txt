##############################################################
# MODULE: AUXILLIARY
##############################################################
# Define directories
set( AUX_DIR    ${SRC}/aux         )
set( AUX_INC    ${AUX_DIR}/include )
set( AUX_SRC    ${AUX_DIR}/src     )

##############################################################
# LIBRARY: AUX
##############################################################
# define name
set( AUX_LIB aux )

# add sources to library
add_library( ${AUX_LIB} INTERFACE  )

# set includes
target_include_directories( ${AUX_LIB}
  INTERFACE $<BUILD_INTERFACE:${AUX_INC}>
  )

install( TARGETS aux DESTINATION ${LIB} )
