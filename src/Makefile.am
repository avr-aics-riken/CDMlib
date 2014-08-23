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

lib_LIBRARIES = libCIO.a

libCIO_a_CXXFLAGS = @MPI_CFLAGS@ @TP_CFLAGS@ -I$(top_builddir)/include

libCIO_a_SOURCES = \
    cio_ActiveSubDomain.C \
    cio_DFI.C \
    cio_FileInfo.C \
    cio_FilePath.C \
    cio_Unit.C \
    cio_TimeSlice.C \
    cio_Domain.C \
    cio_MPI.C \
    cio_Process.C \
    cio_DFI_BOV.C \
    cio_DFI_SPH.C \
    cio_DFI_AVS.C \
    cio_DFI_PLOT3D.C \
    cio_DFI_VTK.C \
    cio_DFI_Read.C \
    cio_DFI_Write.C \
    cio_TextParser.C \
    cio_interp.f90 \
    cio_interp_ijkn.h \
    cio_interp_nijk.h

CIO_includedir = $(includedir)
CIO_include_HEADERS = \
    $(top_builddir)/include/cio_Array.h \
    $(top_builddir)/include/cio_TypeArray.h \
    $(top_builddir)/include/cio_ActiveSubDomain.h \
    $(top_builddir)/include/cio_DFI.h \
    $(top_builddir)/include/cio_FileInfo.h \
    $(top_builddir)/include/cio_FilePath.h \
    $(top_builddir)/include/cio_Unit.h \
    $(top_builddir)/include/cio_TimeSlice.h \
    $(top_builddir)/include/cio_Domain.h \
    $(top_builddir)/include/cio_MPI.h \
    $(top_builddir)/include/cio_Process.h \
    $(top_builddir)/include/cio_DFI_BOV.h \
    $(top_builddir)/include/cio_DFI_SPH.h \
    $(top_builddir)/include/cio_DFI_AVS.h \
    $(top_builddir)/include/cio_DFI_PLOT3D.h \
    $(top_builddir)/include/cio_DFI_VTK.h \
    $(top_builddir)/include/cio_Define.h \
    $(top_builddir)/include/cio_Version.h \
    $(top_builddir)/include/cio_TextParser.h \
    $(top_builddir)/include/cio_PathUtil.h \
    $(top_builddir)/include/cio_endianUtil.h

CIO_inline_includedir = $(includedir)/inline
CIO_inline_include_HEADERS = \
    $(top_builddir)/include/inline/cio_Array_inline.h \
    $(top_builddir)/include/inline/cio_DFI_inline.h \
    $(top_builddir)/include/inline/cio_Plot3d_inline.h


