###################################################################################
#
# FFV-C : Frontflow / violet Cartesian
#
# Copyright (c) 2007-2011 VCAD System Research Program, RIKEN. 
# All rights reserved.
#
# Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo. 
# All rights reserved.
#
# Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved. 
#
###################################################################################

bin_PROGRAMS = fconv 

fconv_CXXFLAGS = \
  -I$(top_builddir)/include \
  -I$(top_builddir)/tools/fconv/include \
  @MPI_CFLAGS@ \
  @CPM_CFLAGS@ \
  @TP_CFLAGS@ \
  @CDM_BUFFER_TUNING@
  
fconv_SOURCES = \
  main.C \
  InputParam.C \
  conv.C \
  convMx1.C \
  convMxM.C \
  convMxN.C \
  convOutput.C \
  convOutput_SPH.C \
  convOutput_BOV.C \
  convOutput_AVS.C \
  convOutput_VTK.C \
  convOutput_PLOT3D.C

fconv_includedir = $(includedir)
fconv_include_HEADERS = \
  $(top_builddir)/tools/fconv/include/conv_Define.h \
  $(top_builddir)/tools/fconv/include/InputParam.h \
  $(top_builddir)/tools/fconv/include/conv.h \
  $(top_builddir)/tools/fconv/include/convMx1.h \
  $(top_builddir)/tools/fconv/include/convMxM.h \
  $(top_builddir)/tools/fconv/include/convMxN.h \
  $(top_builddir)/tools/fconv/include/convOutput.h \
  $(top_builddir)/tools/fconv/include/convOutput_SPH.h \
  $(top_builddir)/tools/fconv/include/convOutput_BOV.h \
  $(top_builddir)/tools/fconv/include/convOutput_AVS.h \
  $(top_builddir)/tools/fconv/include/convOutput_VTK.h \
  $(top_builddir)/tools/fconv/include/convOutput_PLOT3D.h


fconv_inline_includedir = $(includedir)/inline
fconv_inline_include_HEADERS = \
  $(top_builddir)/tools/fconv/include/inline/conv_inline.h \
  $(top_builddir)/tools/fconv/include/inline/convMx1_inline.h \
  $(top_builddir)/tools/fconv/include/inline/conv_plot3d_inline.h

fconv_LDADD = \
  -L$(top_builddir)/src \
  @CPM_LDFLAGS@ \
  @CDM_LIBS@ \
  @TP_LDFLAGS@ \
  @MPI_LDFLAGS@ \
  @MPI_LIBS@ 
  

EXTRA_DIST = Makefile_hand depend.inc 
