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

#                                               -*- cmake -*-
# Process this file with cmake to produce a configure script.

cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING "")
set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "" CACHE STRING "")
set(CMAKE_Fortran_FLAGS "" CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "/usr/local/CDMlib" CACHE PATH "")
set(CDM_BUFFER_TUNING " " CACHE STRING "")

project(CDMlib)

set(PACKAGE "CDMlib")
set(PACKAGE_BUGREPORT "keno@riken.jp")
set(CDMLIB_MAJOR_VERSION 0)
set(CDMLIB_MINOR_VERSION 1)
set(CDMLIB_PATCH_VERSION 0)
set(PACKAGE_VERSION
${CDMLIB_MAJOR_VERSION}.${CDMLIB_MINOR_VERSION}.${CDMLIB_PATCH_VERSION})
set(PACKAGE_NAME ${PACKAGE})
set(PACKAGE_STRING "${PACKAGE} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME ${PACKAGE})
set(VERSION ${PACKAGE_VERSION})

# revision No.
set (CDM_REVISION "20140823_1700")

#
# configure include/cdm_Version.h
#
configure_file("${PROJECT_SOURCE_DIR}/include/cdm_Version.h.in"
               "${PROJECT_SOURCE_DIR}/include/cdm_Version.h" @ONLY)

#
# create config.h from config.h.in
#
include(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
include(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
include(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
include(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
include(${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
include(${CMAKE_ROOT}/Modules/CheckSymbolExists.cmake)

#check_type_size("int" SIZEOF_INT)
#message("int = ${SIZEOF_INT}")
check_function_exists("gethostname" HAVE_GETHOSTNAME)
check_include_file("inttypes.h" HAVE_INTTYPES_H)
check_include_file("memory.h" HAVE_MEMORY_H)
check_function_exists("memset" HAVE_MEMSET)
check_function_exists("mkdir" HAVE_MKDIR)
#check_include_file("stdbool.h" HAVE_STDBOOL_H)
#if (${HAVE_STDBOOL_H} MATCHES "1")
#endif (${HAVE_STDBOOL_H} MATCHES "1")
check_include_file("stdint.h" HAVE_STDINT_H)
check_include_file("stdlib.h" HAVE_STDLIB_H)
check_function_exists("strcasecmp" HAVE_STRCASECMP)
check_function_exists("strerror" HAVE_STRERROR)
check_include_file("strings.h" HAVE_STRINGS_H)
check_include_file("string.h" HAVE_STRING_H)
check_include_file("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
check_include_file("unistd.h" HAVE_UNISTD_H)
check_type_size("_Bool" SIZEOF__BOOL)
#message("_Bool = ${SIZEOF__BOOL}") 
if( ${SIZEOF__BOOL} MATCHES 1)
  set(HAVE__BOOL 1)
endif( ${SIZEOF__BOOL} MATCHES 1)

check_include_file("stdlib.h" HAVE_STDC_HEADERS1)
check_include_file("stddef.h" HAVE_STDC_HEADERS2)
#message("HAVE_STDC_HEADERS1 = ${HAVE_STDC_HEADERS1}")
#message("HAVE_STDC_HEADERS2 = ${HAVE_STDC_HEADERS2}")
if(${HAVE_STDC_HEADERS1} MATCHES 1)
if(${HAVE_STDC_HEADERS2} MATCHES 1)
    set(STDC_HEADERS 1)
endif(${HAVE_STDC_HEADERS2} MATCHES 1)
endif(${HAVE_STDC_HEADERS1} MATCHES 1)

configure_file("${PROJECT_SOURCE_DIR}/config.h.cmake"
               "${PROJECT_BINARY_DIR}/config.h")


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)
#message("${CMAKE_CXX_FLAGS}")
#string(REPLACE "-O3" " " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS})

#
# default prefix setting
#
#set(PREFIX "/usr/local/CDMlib" CACHE PATH "Installation directory for
#libraries")
#message("PREFIX = ${PREFIX}")
#set(CMAKE_INSTALL_PREFIX ${PREFIX})
message("CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

#
# C++
#

#set(CXX " " CACHE PATH "C++ compiler")
#set(CXX_FLAGS "" CACHE STRING "C++ compiler flags")

#if (NOT ${CXX} MATCHES " ")
#  set(CMAKE_CXX_COMPILER ${CXX})
#endif (NOT ${CXX} MATCHES " ")

message("CAMKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

if (${CMAKE_CXX_FLAGS} MATCHES "-")
  set(CMAKE_CXX_FLAGS_RELEASE ${CAME_CXX_FLAGS} )
  set(CMAKE_CXX_FLAGS_RELEASE 
      "${CMAKE_CXX_FLAGS_RELEASE} -DHAVE_CONFIG_H")
else() 
  if(NOT ${CMAKE_CXX_FLAGS} MATCHES "-")
     set(CMAKE_CXX_FLAGS_RELEASE "-g -O2 -DNDEBUG -DHAVE_CONFIG_H")
  endif(NOT ${CMAKE_CXX_FLAGS} MATCHES "-")
endif (${CMAKE_CXX_FLAGS} MATCHES "-")

message("CXX = ${CMAKE_CXX_COMPILER}")
message("CAMKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("CXX_FLAGS = ${CXX_FLAGS}")
#if (NOT DEFINED ${CXX_FLAGS})
#  message("aiueo")
#endif (NOT DEFINED ${CXX_FLAGS})


#
# FORTRAN
#
enable_language(Fortran OPTIONAL)
#set(F90 " " CACHE PATH "Fortran compiler")
#set(F90_FLAGS " " CACHE PATH "Fortran compiler flags")

#if (NOT ${F90} MATCHES " ")
#  set(CMAKE_Fortran_COMPILER ${F90}) 
#endif (NOT ${F90} MATCHES " ")
#message(${CMAKE_Fortran_COMPILER})

if (${CMAKE_Fortran_FLAGS} MATCHES "-")
  set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DHAVE_CONFIG_H")
else()
  if (NOT ${CMAKE_Fortran_FLAGS} MATCHES "-")
    set(CMAKE_Fortran_FLAGS_RELEASE "-g -O2 -DNDEBUG -DHAVE_CONFIG_H")
  endif (NOT ${CMAKE_Fortran_FLAGS} MATCHES "-")
endif (${CMAKE_Fortran_FLAGS} MATCHES "-")

#if (CMAKE_Fortran_COMPILER MATCHES "ifort")
#  # ifort
#   message("use intel fortran!")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#   set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
#endif (CMAKE_Fortran_COMPILER MATCHES "ifort")

#
# MPI environment
#
set(CDM_CXX ${CMAKE_CXX_COMPILER})
message("CDM_CXX = ${CDM_CXX}")
if (${CDM_CXX} MATCHES "mpi")
  set(enable_mpi "wrappermpi")
endif (${CDM_CXX} MATCHES "mpi")
message("enable_mpi= ${enable_mpi}")
#
# MPI library
#
option(with_MPI "use MPI library" ON)
set(with_ompi " " CACHE PATH "Specify OpenMPI install directory")

if(NOT with_MPI MATCHES ON)
  set(CXX_FLAGS "-D_CDM_WITHOUT_MPI_")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_FLAGS}")
else ()
  if(NOT ${with_ompi} MATCHES " ")
    set(MPI_DIR ${with_ompi})
    set(MPI_CFLAGS "${MPI_DIR}/include")
    set(MPI_LDFLAGS "${MPI_DIR}/lib")
    set(MPI_LIBS "mpi")
    set(enable_mpi "ompi")
    set(MPI_LIBS2 "")
  endif(NOT ${with_ompi} MATCHES " ")
endif(NOT with_MPI MATCHES ON)
message("CXX_FLAGS = ${CMAKE_CXX_FLAGS_RELEASE}")
message("MPI_DIR = ${MPI_DIR}")
message("MPI_LIBS = ${MPI_LIBS}")
message("enable_mpi = ${enable_mpi}")

#
# Warning
#
if(with_MPI MATCHES ON)
  if(${enable_mpi} MATCHES "wrappermpi")
    if(${with_ompi} MATCHES " ")
       message("  Warning: OpenMPI library is not specified.")
       message("           If compilation error occurs about MPI, specify
-Dwith_ompi option")
    endif(${with_ompi} MATCHES " ")
  endif(${enable_mpi} MATCHES "wrappermpi")
endif(with_MPI MATCHES ON)

#
# CPMlib
#
set(CPM_DIR " " CACHE PATH "Specify CPMlib installed directory")
set(CPM_LIBS " " CACHE STRING "")
message("CPM_DIR = ${CPM_DIR}")
if (${CPM_DIR} MATCHES " ")
#  message( "Error: TextParser library is not specified.")
#  message( "       Specify -DCPM_DIR option")
#  message( SEND_ERROR "")

else()

set(CPM_LIBS "CPM")
execute_process(
  COMMAND ${CPM_DIR}/bin/cpm-config --cflags
  OUTPUT_VARIABLE CPM_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
#string(REPLACE "-I" "" CPM_CFLAGS ${CPM_CFLAGS})
set(CPM_CFLAGS "${CPM_DIR}/include")
#if(NOT ${CDM_BUFFER_TUNING} MATCHES " ")
#  set(CPM_CFLAGS "${CPM_CFLAGS} -CDM_BUFFER_TUNING=${CDM_BUFFER_TUNING}")
#endif(NOT ${CDM_BUFFER_TUNING} MATCHES " ")
message(${CPM_CFLAGS})

execute_process(
  COMMAND ${CPM_DIR}/bin/cpm-config --libs
  OUTPUT_VARIABLE CPM_LDFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
#string(REPLACE "-L" "" CPM_LDFLAGS ${CPM_LDFLAGS})
message(${CPM_LDFLAGS})
if(CPM_LDFLAGS MATCHES "mpi_cxx")
  set(MPI_LIBS2 "mpi_cxx")
endif(CPM_LDFLAGS MATCHES "mpi_cxx")
set(CPM_LDFLAGS "${CPM_DIR}/lib")

endif (${CPM_DIR} MATCHES " ")

#
# TextParser environment
#
set(TP_DIR " " CACHE PATH "Specify TextParser install directory")
message("TP_DIR = ${TP_DIR}")

if(${TP_DIR} MATCHES " ")
  message("  Error: TextParser library is not specified.")
  message("         Specify -DTP_DIR option.")
  message( SEND_ERROR "")
endif(${TP_DIR} MATCHES " ")

execute_process(
  COMMAND ${TP_DIR}/bin/tp-config --cflags
  OUTPUT_VARIABLE TP_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REPLACE "-I" "" TP_CFLAGS ${TP_CFLAGS})

#if( NOT ${CDM_BUFFER_TUNING} MATCHES " ")
#  set(TP_CFLAGS "${TP_CFLAGS} -CDM_BUFFER_TUNING=${CDM_BUFFER_TUNING}")
#endif( NOT ${CDM_BUFFER_TUNING} MATCHES " ")



message("TP_CFLAGS = ${TP_CFLAGS}")

execute_process(
  COMMAND ${TP_DIR}/bin/tp-config --libs
#  COMMAND grep '-lTP'
  OUTPUT_VARIABLE TP_LDFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
#string(REPLACE "-L" "" TP_LDFLAGS ${TP_LDFLAGS})

message("TP_LDFLAGS = ${TP_LDFLAGS}")
set(TP_LIBS "" CACHE STRING "")
if(${TP_LDFLAGS} MATCHES "-lTPmpi")
  set (TP_LIBS "TPmpi")
  message("TP_LIBS = ${TP_LIBS}")
else()
  if(${TP_LDFLAGS} MATCHES "-lTP")
    set(TP_LIBS "TP")
    message("TP_LIBS = ${TP_LIBS}")
  endif(${TP_LDFLAGS} MATCHES "-lTP")
endif(${TP_LDFLAGS} MATCHES "-lTPmpi")
set(TP_LDFLAGS "${TP_DIR}/lib")
 

#
# frm tools
#
option(with_frm "install frm tool" OFF)
if(${CMAKE_CROSSCOMPILING})
  message("  Warning : frm tool is not build due to cross-compilation.")
endif(${CMAKE_CROSSCOMPILING})

#
# CDM special flags
#

message(${CMAKE_INSTALL_PREFIX})
set(CDM_CFLAGS "${CMAKE_INSTALL_PREFIX}/include")
set(CDM_LDFLAGS "${CMAKE_INSTALL_PREFIX}/lib")
set(CDM_LIBS "CDM")

if( NOT ${CDM_BUFFER_TUNING} MATCHES " ")
  set(BUFFER_TUNING_FLAGS "-DCDM_BUFFER_MB_SIZE=${CDM_BUFFER_TUNING}")
  set(CDM_CFLAGS "${CDM_CFLAGS} ${BUFFER_TUNING_FLAGS}")
else()
  set(BUFFER_TUNING_FLAGS " ")
endif( NOT ${CDM_BUFFER_TUNING} MATCHES " ")

message(${CDM_CFLAGS})
message(${CDM_LDFLAGS})
message(${CDM_LIBS})


#
# MAKE_SUB_DIRS
#
make_directory(src)
make_directory(doc)

# install CDMlib
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/doc)

# install frm
message("with_frm = ${with_frm}")
if(${with_frm} MATCHES ON)
#  configure_file("${PROJECT_SOURCE_DIR}/tools/frm/src/CMakeLists.txt.cmake"
#                 "${PROJECT_SOURCE_DIR}/tools/frm/src/CMakeLists.txt" @ONLY)
  make_directory("${PROJECT_SOURCE_DIR}/tools/frm/src")
  add_subdirectory("${PROJECT_SOURCE_DIR}/tools/frm/src")
endif(${with_frm} MATCHES ON)

# install fconv
message("CPM_DIR = ${CPM_DIR}")
if(${with_MPI} MATCHES ON)
  if( NOT ${CPM_DIR} MATCHES " ")
     make_directory("${PROJECT_SOURCE_DIR}/tools/fconv/src")
     add_subdirectory("${PROJECT_SOURCE_DIR}/tools/fconv/src")
  endif( NOT ${CPM_DIR} MATCHES " ")
endif(${with_MPI} MATCHES ON)


#
# configure include/cdm_Version.h
#
set(CMAKE_EXEC_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CMAKE_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)

set(CMAKE_CDM_CFLAGS "-I${CDM_CFLAGS}")
if(NOT ${with_ompi} MATCHES " ")
  set(CMAKE_MPI_CFLAGS "-I${MPI_CFLAGS}")
  set(CMAKE_MPI_LDFLAGS "-L${MPI_LDFLAGS}")
  set(CMAKE_MPI_LIBS "-l${MPI_LIBS}")
else()
  set(CMAKE_MPI_CFLAGS "")
  set(CMAKE_MPI_LDFLAGS "")
  set(CMAKE_MPI_LIBS "")
endif(NOT ${with_ompi} MATCHES " ")
set(CMAKE_TP_CFLAGS "-I${TP_CFLAGS}")
set(CMAKE_CDM_LDFLAGS "-L${CDM_LDFLAGS}") 

set(CMAKE_CDM_LIBS "-l${CDM_LIBS}")
set(CMAKE_TP_LDFLAGS "-L${TP_LDFLAGS}")
set(CMAKE_TP_LIBS "-l${TP_LIBS}")
configure_file("${PROJECT_SOURCE_DIR}/cdm-config.cmake"
               "${PROJECT_SOURCE_DIR}/cdm-config" @ONLY)

execute_process(
  COMMAND chmod a+x cdm-config
)

install(FILES "${PROJECT_SOURCE_DIR}/cdm-config" 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE  WORLD_READ WORLD_EXECUTE
)

install(FILES ${PROJECT_SOURCE_DIR}/AUTHORS 
              ${PROJECT_SOURCE_DIR}/ChangeLog 
              ${PROJECT_SOURCE_DIR}/LICENSE 
              ${PROJECT_SOURCE_DIR}/COPYING 
              ${PROJECT_SOURCE_DIR}/README.md 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share
)

add_custom_target(distclean
   COMMAND /bin/bash ./clean.sh
)


#uninstall target
add_custom_target(uninstall
COMMAND rm -rf ${CMAKE_INSTALL_PREFIX}/bin
               ${CMAKE_INSTALL_PREFIX}/doc
               ${CMAKE_INSTALL_PREFIX}/include
               ${CMAKE_INSTALL_PREFIX}/lib
               ${CMAKE_INSTALL_PREFIX}/share
)


#set(distdir ${PACKAGE}-${PACKAGE_VERSION}) 
#set(DIST_ARCHIVES ${distdir}.tar.gz)
#message("distdir = ${distdir}")
#message("DIST_ARCHIVES = ${DIST_ARCHIVES}")
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE
"${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${CDMLIB_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${CDMLIB_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${CDMLIB_PATCH_VERSION}")
include(CPack)
