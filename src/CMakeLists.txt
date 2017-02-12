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
      ${TP_LIB}
      ${NetCDF_LIB}
)


set(cdm_files
    cdm_ActiveSubDomain.C
    cdm_DFI_AVS.C
    cdm_DFI_BOV.C
    cdm_DFI_FUB.C
    cdm_DFI_NETCDF.C
    cdm_DFI_PLOT3D.C
    cdm_DFI_Read.C
    cdm_DFI_SPH.C
    cdm_DFI_VTK.C
    cdm_DFI_Write.C
    cdm_DFI.C
    cdm_Domain.C
    cdm_FieldFileNameFormat.C
    cdm_FieldFileNameFormatElem.C
    cdm_FileInfo.C
    cdm_FilePath.C
    cdm_interp.f90
    cdm_MPI.C
    cdm_Process.C
    cdm_TextParser.C
    cdm_TimeSlice.C
    cdm_Unit.C
    cdm_VisIt.C
)


if(with_MPI)
  add_library(CDMmpi STATIC ${cdm_files})
  target_link_libraries(CDMmpi -lCPMmpi -lTPmpi)
  install(TARGETS CDMmpi DESTINATION lib)
else()
  add_library(CDM STATIC ${cdm_files})
  target_link_libraries(CDM -lCPM -lTP)
  install(TARGETS CDM DESTINATION lib)
endif()


install(FILES
        ${PROJECT_SOURCE_DIR}/include/cdm_ActiveSubDomain.h
        ${PROJECT_SOURCE_DIR}/include/cdm_Array.h
        ${PROJECT_SOURCE_DIR}/include/cdm_Define.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_AVS.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_BOV.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_FUB.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_NETCDF.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_PLOT3D.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_SPH.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI_VTK.h
        ${PROJECT_SOURCE_DIR}/include/cdm_DFI.h
        ${PROJECT_SOURCE_DIR}/include/cdm_Domain.h
        ${PROJECT_SOURCE_DIR}/include/cdm_endianUtil.h
        ${PROJECT_SOURCE_DIR}/include/cdm_FieldFileNameFormat.h
        ${PROJECT_SOURCE_DIR}/include/cdm_FieldFileNameFormatElem.h
        ${PROJECT_SOURCE_DIR}/include/cdm_FILE.h
        ${PROJECT_SOURCE_DIR}/include/cdm_FileInfo.h
        ${PROJECT_SOURCE_DIR}/include/cdm_FilePath.h
        ${PROJECT_SOURCE_DIR}/include/cdm_MPI.h
        ${PROJECT_SOURCE_DIR}/include/cdm_NonUniformDomain.h
        ${PROJECT_SOURCE_DIR}/include/cdm_PathUtil.h
        ${PROJECT_SOURCE_DIR}/include/cdm_Process.h
        ${PROJECT_SOURCE_DIR}/include/cdm_TextParser.h
        ${PROJECT_SOURCE_DIR}/include/cdm_TimeSlice.h
        ${PROJECT_SOURCE_DIR}/include/cdm_TypeArray.h
        ${PROJECT_SOURCE_DIR}/include/cdm_Unit.h
        ${PROJECT_SOURCE_DIR}/include/cdm_VisIt.h
        ${PROJECT_SOURCE_DIR}/include/mpi_stubs.h
        ${PROJECT_BINARY_DIR}/include/cdm_Version.h
        DESTINATION include
)

install(FILES
        ${PROJECT_SOURCE_DIR}/include/inline/cdm_Array_inline.h
        ${PROJECT_SOURCE_DIR}/include/inline/cdm_DFI_inline.h
        ${PROJECT_SOURCE_DIR}/include/inline/cdm_NonUniformDomain_inline.h
        ${PROJECT_SOURCE_DIR}/include/inline/cdm_Plot3d_inline.h
        DESTINATION include/inline
)
