###
#########################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
#########################################################################
###

lib_LIBRARIES = libCDM.a

libCDM_a_CXXFLAGS = @MPI_CFLAGS@ @TP_CFLAGS@ -I$(top_builddir)/include

libCDM_a_SOURCES = \
    cdm_ActiveSubDomain.C \
    cdm_DFI.C \
    cdm_FileInfo.C \
    cdm_FilePath.C \
    cdm_VisIt.C \
    cdm_Unit.C \
    cdm_TimeSlice.C \
    cdm_Domain.C \
    cdm_MPI.C \
    cdm_Process.C \
    cdm_DFI_BOV.C \
    cdm_DFI_SPH.C \
    cdm_DFI_AVS.C \
    cdm_DFI_PLOT3D.C \
    cdm_DFI_VTK.C \
    cdm_DFI_Read.C \
    cdm_DFI_Write.C \
    cdm_TextParser.C \
    cdm_interp.f90 \
    cdm_interp_ijkn.h \
    cdm_interp_nijk.h

CDM_includedir = $(includedir)
CDM_include_HEADERS = \
    $(top_builddir)/include/cdm_Array.h \
    $(top_builddir)/include/cdm_TypeArray.h \
    $(top_builddir)/include/cdm_ActiveSubDomain.h \
    $(top_builddir)/include/cdm_DFI.h \
    $(top_builddir)/include/cdm_FileInfo.h \
    $(top_builddir)/include/cdm_FilePath.h \
    $(top_builddir)/include/cdm_VisIt.h \
    $(top_builddir)/include/cdm_Unit.h \
    $(top_builddir)/include/cdm_TimeSlice.h \
    $(top_builddir)/include/cdm_Domain.h \
    $(top_builddir)/include/cdm_NonUniformDomain.h \
    $(top_builddir)/include/cdm_MPI.h \
    $(top_builddir)/include/cdm_Process.h \
    $(top_builddir)/include/cdm_DFI_BOV.h \
    $(top_builddir)/include/cdm_DFI_SPH.h \
    $(top_builddir)/include/cdm_DFI_AVS.h \
    $(top_builddir)/include/cdm_DFI_PLOT3D.h \
    $(top_builddir)/include/cdm_DFI_VTK.h \
    $(top_builddir)/include/cdm_Define.h \
    $(top_builddir)/include/cdm_Version.h \
    $(top_builddir)/include/cdm_TextParser.h \
    $(top_builddir)/include/cdm_PathUtil.h \
    $(top_builddir)/include/cdm_endianUtil.h

CDM_inline_includedir = $(includedir)/inline
CDM_inline_include_HEADERS = \
    $(top_builddir)/include/inline/cdm_Array_inline.h \
    $(top_builddir)/include/inline/cdm_DFI_inline.h \
    $(top_builddir)/include/inline/cdm_Plot3d_inline.h \
    $(top_builddir)/include/inline/cdm_BOV_inline.h \
    $(top_builddir)/include/inline/cdm_NonUniformDomain_inline.h


