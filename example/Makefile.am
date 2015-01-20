#########################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
#########################################################################

noinst_PROGRAMS = input output output_nuc


DISTCLEANFILES=*~ input output output_nuc
CLEANFILES=

EXTRA_DIST= sample_input.C sample_output_NonUni.C sample_output_Uni.C coord.crd coord_Ascii.crd


input_SOURCES  =  sample_input.C
input_CXXFLAGS = -I$(top_builddir)/include @MPI_CFLAGS@ @TP_CFLAGS@

input_LDADD = \
    -L$(top_builddir)/src -lCDM \
    @MPI_LIBS@ \
    @MPI_LDFLAGS@ \
    @TP_LDFLAGS@


output_SOURCES  = sample_output_Uni.C
output_CXXFLAGS = -I$(top_builddir)/include @MPI_CFLAGS@ @TP_CFLAGS@

output_LDADD = \
    -L$(top_builddir)/src -lCDM \
    @MPI_LIBS@ \
    @MPI_LDFLAGS@ \
    @TP_LDFLAGS@


output_nuc_SOURCES  = sample_output_NonUni.C
output_nuc_CXXFLAGS = -I$(top_builddir)/include @MPI_CFLAGS@ @TP_CFLAGS@

output_nuc_LDADD = \
  -L$(top_builddir)/src -lCDM \
  @MPI_LIBS@ \
  @MPI_LDFLAGS@ \
  @TP_LDFLAGS@


dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

