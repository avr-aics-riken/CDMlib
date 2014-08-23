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

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS = cdm-config
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = cdm-config.in LICENSE include/mpi_stubs.h
dist_data_DATA = AUTHORS COPYING ChangeLog README.md LICENSE
dist_include_HEADERS = include/mpi_stubs.h

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

