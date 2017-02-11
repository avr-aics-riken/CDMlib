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

# - Try to find CDMlib
# Once done, this will define
#
#  CDM_FOUND - system has CDMlib
#  CDM_INCLUDE_DIRS - CDMlib include directories
#  CDM_LIBRARIES - link these to use CDMlib

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(CDM_PKGCONF CDM)

if(CMAKE_PREFIX_PATH)
  set(CDM_CANDIDATE_PATH ${CMAKE_PREFIX_PATH})
  file(GLOB tmp "${CMAKE_PREFIX_PATH}/[Jj][Hh][Pp][Cc][Nn][Dd][Ff]*/")
  list(APPEND CDM_CANDIDATE_PATH ${tmp})
endif()

# Include dir
find_path(CDM_INCLUDE_DIR
  NAMES cdm_Define.h
  PATHS ${CDM_ROOT} ${CDM_PKGCONF_INCLUDE_DIRS} ${CDM_CANDIDATE_PATH}
  PATH_SUFFIXES include
)

# Finally the library itself
find_library(CDM_LIBRARY
  NAMES CDM
  PATHS ${CDM_ROOT} ${CDM_PKGCONF_LIBRARY_DIRS} ${CDM_CANDIDATE_PATH}
  PATH_SUFFIXES lib
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(CDM_PROCESS_INCLUDES CDM_INCLUDE_DIR)
set(CDM_PROCESS_LIBS CDM_LIBRARY)
libfind_process(CDM)
