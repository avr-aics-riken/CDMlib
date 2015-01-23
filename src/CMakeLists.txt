###
#########################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
#########################################################################
###

message("entered ${PROJECT_SOURCE_DIR}/src")

#set(lib_CDM_a_CXXFLAGS "${MPI_CFLAGS} ${TP_CFLAGS}")
#set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_INIT "")
#set(CMAKE_CXX_FLAGS "")
#set(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING "" FORCE)
#set(CMAKE_CXX_FLAGS "-DNDEBUG" CACHE STRING "" FORCE)

#add_definitions(${CMAKE_CXX_FLAGS} "-g -O2")

#configure_file(
#  "${PROJECT_SOURCE_DIR}/config.h.in"
#  "${PROJECT_SOURCE_DIR}/config.h"
#)

#"-Dinclude_name="hoge -E""

include_directories( 
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/src 
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/include/inline
  ${MPI_DIR}/include
  ${TP_CFLAGS}
)

#link_directories(
#  "${MPI_DIR}/lib -lmpi"
#)


set(libCDM_a_SOURCES
    cdm_ActiveSubDomain.C
    cdm_DFI.C
    cdm_FileInfo.C
    cdm_FilePath.C
    cdm_Unit.C
    cdm_TimeSlice.C
    cdm_Domain.C
#    cdm_NonUniformDomain.C
    cdm_MPI.C
    cdm_Process.C
    cdm_DFI_BOV.C
    cdm_DFI_SPH.C
    cdm_DFI_AVS.C
    cdm_DFI_PLOT3D.C
    cdm_DFI_VTK.C
    cdm_DFI_Read.C
    cdm_DFI_Write.C
    cdm_TextParser.C
    cdm_interp.f90
    cdm_VisIt.C
)

add_definitions(${BUFFER_TUNING_FLAGS})

add_library(CDM STATIC ${libCDM_a_SOURCES} ${PROJECT_SOURCE_DIR}/config.h)

set_target_properties(CDM PROPERTIES PUBLIC_HEADER ${PROJECT_SOURCE_DIR}/config.h)

install(TARGETS CDM ARCHIVE DESTINATION  ${CMAKE_INSTALL_PREFIX}/lib)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include DESTINATION
${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.h")

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/inline DESTINATION
${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")
