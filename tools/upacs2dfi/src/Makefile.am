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
# Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved. 
#
###################################################################################

bin_PROGRAMS = upacs2dfi 

upacs2dfi_CXXFLAGS = \
  -I$(top_srcdir)/include \
  -I$(top_srcdir)/tools/upacs2dfi/src \
  @MPI_CFLAGS@ \
  @TP_CFLAGS@ \
  @CDM_BUFFER_TUNING@ 
  
upacs2dfi_SOURCES = \
  fub_DATA.C \
  main.C

upacs2dfi_includedir = $(includedir)
upacs2dfi_include_HEADERS = \
  $(top_srcdir)/tools/upacs2dfi/src/fub_Define.h \
  $(top_srcdir)/tools/upacs2dfi/src/Coord.h \
  $(top_srcdir)/tools/upacs2dfi/src/fub_DATA.h

upacs2dfi_LDADD = \
  -L$(top_builddir)/src \
  @CDM_LIBS@ \
  @TP_LDFLAGS@ \
  @MPI_LDFLAGS@ \
  @MPI_LIBS@ 
  
