#########################################################################
#
# CDMlib netcdf2dfi
# 
# Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
# #######################################################################

bin_PROGRAMS = netcdf2dfi

netcdf2dfi_CXXFLAGS = \
  -I$(top_srcdir)/include \
  -I$(top_srcdir)/tools/netcdf2dfi/include \
  @TP_CFLAGS@ \
  @CDM_BUFFER_TUNING@ \
  @NC_CFLAGS@

netcdf2dfi_SOURCES = \
  main.C \
  netcdf2dfi.C

netcdf2dfi_includedir = $(includedir)
netcdf2dfi_include_HEADERS = \
  $(top_srcdir)/tools/netcdf2dfi/include/netcdf2dfi.h

netcdf2dfi_LDADD = \
  -L$(top_builddir)/src \
  @CDM_LIBS@ \
  @TP_LDFLAGS@ \
  @NC_LDFLAGS@

#EXTRA_DIST = Makefile_hand
