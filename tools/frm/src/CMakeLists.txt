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
       ${PROJECT_SOURCE_DIR}/tools/frm/include
       ${TP_INC}
       ${NetCDF_INC}
       ${HDF_INC}
)


link_directories(
      ${PROJECT_BINARY_DIR}/src
      ${TP_LIB}
      ${NetCDF_LIB}
      ${HDF_LIB}
)


set(frm_files
    ActiveSubDomain.C
    main.C
    Staging.C
)


if(with_MPI)
  add_executable(frm ${frm_files})
  target_link_libraries(frm -lCDMmpi -lTPmpi ${LIB_NETCDF})
  install(TARGETS frm DESTINATION bin)
else()
  add_executable(frm ${frm_files})
  target_link_libraries(frm -lCDM -lTP ${LIB_NETCDF})
  install(TARGETS frm DESTINATION bin)
endif()
