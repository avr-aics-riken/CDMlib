#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
#! /bin/sh

prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@CMAKE_EXEC_PREFIX@
includedir=@CMAKE_INCLUDE_DIR@
libdir=@CMAKE_LIB_DIR@

usage()
{
    cat <<EOF

Usage: cdm-config [OPTION]

Known values for OPTION are:

  --cxx       print C++ compiler command
  --cflags    print C/C++ pre-processor and compiler flags
  --libs      print library linking information for C++ program
  --help      display this help and exit
  --version   output version information

EOF

    exit $1
}

if test $# -eq 0; then
    usage 1
fi

cflags=false
libs=false

while test $# -gt 0; do
    case "$1" in
    -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) optarg= ;;
    esac

    case "$1" in
    --version)
cat <<EOF

CDMlib : Cartesian Data Management library  Version : @VERSION@ : @CDM_REVISION@

Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.

EOF
      exit 0
      ;;

    --help)
      usage 0
      ;;

    --cxx)
      echo @CDM_CXX@
      ;;

    --cflags)
      echo @CMAKE_CDM_CFLAGS@ @CMAKE_MPI_CFLAGS@ @CMAKE_TP_CFLAGS@
      ;;

    --libs)
      echo @CMAKE_CDM_LDFLAGS@ @CMAKE_CDM_LIBS@ @CMAKE_MPI_LDFLAGS@ @CMAKE_MPI_LIBS@ @CMAKE_TP_LDFLAGS@ @CMAKE_TP_LIBS@
      ;;

    *)
      usage
      exit 1
      ;;
    esac
    shift
done

exit 0
