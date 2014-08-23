#########################################################################
#
# CIOlib Staging 
# 
# Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
# #######################################################################

bin_PROGRAMS = frm

frm_CXXFLAGS = \
  -I$(top_builddir)/include \
  -I$(top_builddir)/tools/frm/include \
  @TP_CFLAGS@ 

frm_SOURCES = \
  Staging.C \
  ActiveSubDomain.C \
  main.C

frm_includedir = $(includedir)
frm_include_HEADERS = \
  $(top_builddir)/tools/frm/include/ActiveSubDomain.h \
  $(top_builddir)/tools/frm/include/Staging_Define.h \
  $(top_builddir)/tools/frm/include/Staging_Utility.h \
  $(top_builddir)/tools/frm/include/frm_EndianUtil.h \
  $(top_builddir)/tools/frm/include/frm_Version.h

frm_LDADD = \
  -L$(top_builddir)/src \
  @CIO_LIBS@ \
  @TP_LDFLAGS@

EXTRA_DIST = Makefile_hand
