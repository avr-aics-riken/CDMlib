###################################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
#
#-*- CMAKE -*-#
# Process this file with cmake to produce a configure script.

# Following parameters must be given at cmake command line

# Build options, see Readme.md for details.
#
# -D INSTALL_DIR=install_directory
#
# -D with_MPI={yes|no}
#
# -D with_TP=installed_directory
#
# -D with_CPM={installed_directory|no}
#
# -D with_NetCDF=installed_directory
#
# -D with_HDF=installed_directory
#
# -D with_util={yes|no}
#
# -D with_example={no|yes}
#
# -D enable_BUFFER_SIZE=<size MB>
#

cmake_minimum_required(VERSION 2.6)


#######
# Build options
#######

option (with_util "Enable utility" "ON")
option (with_MPI "Enable MPI" "ON")
option (with_example "Compiling examples" "OFF")
#option (with_NetCDF "Enable NetCDF" "OFF")
#option (with_HDF "Enable NetCDF(HDF)" "OFF")
#option (with_CPM "Enable CPMlib" "OFF")
#option (with_TP "Enable TextParser" "OFF")
option (enable_BUFFER_SIZE "Enable Buffer Size option" "OFF")


#######
# Project setting
#######

project(CDMlib C CXX Fortran)

if(INSTALL_DIR STREQUAL "")
  set(INSTALL_DIR "/usr/local/CDMlib")
endif()

# project()でCMAKE_INSTALL_PREFIX=/usr/localとなるので、再定義
set (CMAKE_INSTALL_PREFIX "${INSTALL_DIR}")

#set(CMAKE_BUILD_TYPE Release)



#######
set(PROJECT_VERSION "1.1.3")
set(LIB_REVISION "20170705_2200")
#######



#######
# Cmake module
#######

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileOptionSelector.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/cdm_package.cmake)




#######
# Options
#######

# Fortran free form
FreeForm()

AddOptimizeOption()

#checkOpenMP()

if(enable_BUFFER_SIZE)
  if(enable_BUFFER_SIZE STREQUAL "no" )
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCDM_BUFFER_MB_SIZE=${enable_BUFFER_SIZE}")
  endif()
endif()


# Serial flag
if(NOT with_MPI)
  SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -DDISABLE_MPI")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MPI")
endif()



#######
# Display options
#######

message(" ")
message( STATUS "Destination PATH: "               ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_C_COMPILER: "               ${CMAKE_C_COMPILER})
message( STATUS "CMAKE_C_FLAGS: "                  ${CMAKE_C_FLAGS})
message( STATUS "CMAKE_CXX_COMPILER: "             ${CMAKE_CXX_COMPILER})
message( STATUS "CMAKE_CXX_FLAGS: "                ${CMAKE_CXX_FLAGS})

if(enable_Fortran)
  message( STATUS "CMAKE_Fortran_COMPILER: "         ${CMAKE_Fortran_COMPILER})
  message( STATUS "CMAKE_Fortran_FLAGS: "            ${CMAKE_Fortran_FLAGS})
endif()

message(" ")
message( STATUS "Utilities         : "      ${with_util})
message( STATUS "MPI support       : "      ${with_MPI})
message( STATUS "TextParser support: "      ${with_TP})
message( STATUS "CPMlib support    : "      ${with_CPM})
message( STATUS "NetCDF support    : "      ${with_NetCDF})
message( STATUS "HDF support       : "      ${with_HDF})
message( STATUS "Example           : "      ${with_example})
message( STATUS "BUFFER_SIZE       : "      ${enable_BUFFER_SIZE})
message(" ")

if(USE_F_TCS STREQUAL "YES")
else()
  if(with_MPI)
    include(FindMPI)

    if(MPI_CXX_FOUND)
      message("MPI-CXX-CMP  : ${MPI_CXX_COMPILER}")
      message("MPI-CXX-FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
      message("MPI-CXX-INC  : ${MPI_CXX_INCLUDE_PATH}")
      message("MPI-LNK-FLAGS: ${MPI_CXX_LINK_FLAGS}")
      message("MPI-LIBS     : ${MPI_CXX_LIBRARIES}")
    endif(MPI_CXX_FOUND)
  endif()
endif()


message(" ")
message("PROJECT : CDMlib " ${PROJECT_VERSION} " : " ${LIB_REVISION})
message(" ")



#######
# include PATH
#######

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/inline
    )


#######
# Find Libraries to depend on
#######

# TextParser
IF(with_TP)
  SET(TP_DIR "${with_TP}")
  SET(TP_INC "${TP_DIR}/include")
  SET(TP_LIB "${TP_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find TextParser.")
ENDIF()


# CPMlib
IF(with_CPM STREQUAL "no")
ELSE()
  SET(CPM_DIR "${with_CPM}")
  SET(CPM_INC "${CPM_DIR}/include")
  SET(CPM_LIB "${CPM_DIR}/lib")
ENDIF()

# NetCDF
IF(with_NetCDF STREQUAL "no")
else()
  IF(with_HDF)
    SET(HDF_DIR "${with_HDF}")
    SET(HDF_INC "${HDF_DIR}/include")
    SET(HDF_LIB "${HDF_DIR}/lib")
    SET(NetCDF_DIR "${with_NetCDF}")
    SET(NetCDF_INC "${NetCDF_DIR}/include")
    SET(NetCDF_LIB "${NetCDF_DIR}/lib")
    SET(LIB_NETCDF "-lnetcdf -lhdf5_hl -lhdf5 -lz")
  ELSE()
    SET(NetCDF_DIR "${with_NetCDF}")
    SET(NetCDF_INC "${NetCDF_DIR}/include")
    SET(NetCDF_LIB "${NetCDF_DIR}/lib")
    SET(LIB_NETCDF "-lnetcdf")
  ENDIF()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_WITH_NETCDF4_")
ENDIF()


## Check frm tool
if(with_util)
  if(CDM_CROSS_OPTION)
    message(" ")
    message("frm        : not build because cross compilation is assumed")
  else()
    message("frm        : build")
  endif()
endif()

if(with_MPI)
  message("fconv      : build")
  message("upacs2dfi  : build")
else()
  message("fconv      : not build because with_MPI=no")
  message("upacs2dfi  : not build because with_MPI=no")
endif()

if(with_NetCDF STREQUAL "no")
  message("netcdf2dfi : not build because with_NetCDF=no")
else()
  message("netcdf2dfi : build")
endif()

add_definitions(-DHAVE_CONFIG_H)



#######
# Example
#######

if (with_example)
  enable_testing()
  message("Code test will be performed")
else()
  message("Code test is skipped")
endif()



#######
# Check header files
#######

include(CheckIncludeFiles)

CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)



#######
# Sub Directories
#######

add_subdirectory(src)

if(with_util)
  add_subdirectory(tools)
endif()

add_subdirectory(doc)

if(with_example)
  add_subdirectory(example)
endif()


#######
# configure files
#######

configure_file( config.h.in.cmake config.h )
configure_file( include/cdm_Version.h.in ${PROJECT_BINARY_DIR}/include/cdm_Version.h @ONLY)


#######
# Install files
#######

install(FILES  ChangeLog.md
               Readme.md
               License.txt
               ${PROJECT_BINARY_DIR}/CMakeCache.txt
        DESTINATION share)


#######
# Installer
#######

cdm_package()
