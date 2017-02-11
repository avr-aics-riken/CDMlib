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
       ${TP_INC}
       ${NetCDF_INC}
)


link_directories(
      ${PROJECT_BINARY_DIR}/src
      ${TP_LIB}
      ${NetCDF_LIB}
)


set(upacs2dfi_files
    fub_DATA.C
    main.C
)


add_library(upacs2dfi STATIC ${upacs2dfi_files})
target_link_libraries(upacs2dfi -lCDMmpi -lTPmpi)
install(TARGETS upacs2dfi DESTINATION lib)
