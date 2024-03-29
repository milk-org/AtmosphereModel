set(LIBNAME "milkAtmosphereModel")

message("")
message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")

# list source files (.c) other than modulename.c
set(SOURCEFILES
	AirMixture_ria.c
	AtmosphereModel_build_stdAtmModel.c
    AtmosphereModel_Create_from_CONF.c
    AtmosphereModel_H2O_Saturation.c
    AtmosphereModel_RefractionPath.c
    AtmosphereModel_save_stdAtmModel.c
    AtmosphereModel_stdAtmModel_ria.c
    nrlmsise-00.c
    nrlmsise-00_data.c
    nrlmsise-00_test.c
    )

# list include files (.h) that should be installed on system
set(INCLUDEFILES
	AtmosphereModel_Create_from_CONF.h
	)

# list scripts that should be installed on system
set(SCRIPTFILES
	)

set(LINKLIBS
	CLIcore
)


# DEFAULT SETTINGS
# Do not change unless needed
# =====================================================================

# SRCNAME is current directory (last subdir)
set(NAME0 "${CMAKE_CURRENT_SOURCE_DIR}")
string(REPLACE "/" " " NAME1 ${NAME0})
string(REGEX MATCH "[a-zA-Z0-9_]+$" SRCNAME "${NAME0}")
message(" SRCNAME = ${SRCNAME}")

# LIBNAME constructed by removing spaces
#string(REPLACE "_" "" LIBNAME "${PROJECT_NAME}${SRCNAME}")
#message(" LIBNAME = ${LIBNAME}")

project(lib_${LIBNAME}_project)

# Library can be compiled from multiple source files
# Convention: the main souce file is named <libname>.c
#
add_library(${LIBNAME} SHARED ${SRCNAME}.c ${SOURCEFILES})

target_include_directories(${LIBNAME} PUBLIC ${PROJECT_SOURCE_DIR}/.. ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${LIBNAME} PRIVATE ${LINKLIBS})

install(TARGETS ${LIBNAME} DESTINATION lib)
install(FILES ${SRCNAME}.h ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(PROGRAMS ${SCRIPTFILES} DESTINATION bin)
