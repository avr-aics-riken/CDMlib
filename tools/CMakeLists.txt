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

if(CDM_CROSS_OPTION)
else()
  add_subdirectory(frm/src)
endif()

if(with_MPI)
  add_subdirectory(fconv/src)
  add_subdirectory(upacs2dfi/src)
endif()

if(with_NetCDF STREQUAL "OFF")
else()
  add_subdirectory(netcdf2dfi/src)
endif()
