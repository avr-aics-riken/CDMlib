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

# cdm_Version.hはPROJECT_BINARY_DIR
include_directories(
       ${PROJECT_BINARY_DIR}/include
       ${PROJECT_SOURCE_DIR}/include
       ${PROJECT_SOURCE_DIR}/include/inline
       ${PROJECT_SOURCE_DIR}/tools/ntcdf2dfi/include
       ${TP_INC}
       ${NetCDF_INC}
)


link_directories(
      ${PROJECT_BINARY_DIR}/src
      ${TP_LIB}
      ${NetCDF_LIB}
)


set(netcdf2dfi_files
    main.C
    netcdf2dfi.C
)


add_executable(netcdf2dfi ${netcdf2dfi_files})
target_link_libraries(netcdf2dfi -lCDMmpi -lTPmpi)
install(TARGETS netcdf2dfi DESTINATION bin)
