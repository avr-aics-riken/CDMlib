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
       ${PROJECT_SOURCE_DIR}/tools/fconv/include
       ${PROJECT_SOURCE_DIR}/tools/fconv/include/inline
       ${TP_INC}
       ${CPM_INC}
       ${NetCDF_INC}
)


link_directories(
      ${PROJECT_BINARY_DIR}/src
      ${TP_LIB}
      ${CPM_LIB}
      ${NetCDF_LIB}
)


set(fconv_files
    conv.C
    convMx1.C
    convMxM.C
    convMxN.C
    convOutput_AVS.C
    convOutput_BOV.C
    convOutput_FUB.C
    convOutput_NETCDF.C
    convOutput_PLOT3D.C
    convOutput_SPH.C
    convOutput_VTK.C
    convOutput.C
    InputParam.C
    main.C
)


add_library(fconv STATIC ${fconv_files})
target_link_libraries(fconv -lCDMmpi -lCPMmpi -lTPmpi)
install(TARGETS fconv DESTINATION lib)
