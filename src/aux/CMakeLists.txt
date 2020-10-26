##############################################################
# MODULE: AUXILLIARY
##############################################################
# Define directories
set( AUX_DIR    ${SRC}/aux         )
set( AUX_INC    ${AUX_DIR}/include )
set( AUX_SRC    ${AUX_DIR}/src     )

# Define source files
set( AUX_MAIN
  ${AUX_SRC}/bstrlib.c
  )

##############################################################
# LIBRARY: AUX
##############################################################
# define name
set( AUX_LIB aux )

# add sources to library
#add_library( ${AUX_LIB} INTERFACE  )
add_library( ${AUX_LIB} STATIC ${AUX_MAIN} )

# set includes
#target_include_directories( ${AUX_LIB}
#  INTERFACE $<BUILD_INTERFACE:${AUX_INC}>
#  )

target_include_directories( ${AUX_LIB}
  PUBLIC $<BUILD_INTERFACE:${AUX_INC}>
  PRIVATE ${AUX_SRC}
)

install( TARGETS aux DESTINATION ${LIB} )
