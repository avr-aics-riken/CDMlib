# CDMlib - Cartesian Data Management library

* Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN. All rights reserved.
* Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University. All rights reserved.


## OUTLINE

CDMlib provides following functions for reading and writing files of the Cartesian and Non-uniform Cartesian data structures.

  - Management of distributed files
  - Support file format : SPH, BOV, PLOT3D, NetCDF
  - Restart from previous calculated data (standard)
  - Restart from previous calculated data of the number of different process
  - Restart from previous coarse data with interpolation
  - Staging helper
  - File converter utility {SPH, BOV, PLOT3D, NetCDF} >> {VTK, AVS, SPH, BOV, PLOT3D}


## SOFTWARE REQUIREMENT
- Cmake
- MPI library
- TextParser
- CPMlib(option, required to build tools, not for library itself)
- NetCDF4 library(option)
- HDF library(option)


## INGREDIENTS
~~~
ChangeLog         History of development
License.txt       License to apply
Readme.md         This document, including the description of build
cmake/            Modules of cmake
doc/              Document
example/          Example source codes
include/          Header files
src/              Source codes
tools/frm         File Rank Mapper
tools/fconv       File Converter
tools/netcdf2dfi  DFI generator for NetCDF
tools/upacs2dfi   DFI generator for UPACS code
~~~


## HOW TO BUILD

### Build

~~~
$ export CDM_HOME=/hogehoge
$ mkdir build
$ cd build
$ cmake [options] ..
$ make
$ sudo make install
~~~


### Options

`-D INSTALL_DIR=` *Install_directory*

>  Specify the directory that this library will be installed. Built library is installed at `install_directory/lib` and the header files are placed at `install_directory/include`. The default install directory is `/usr/local/CDMlib`.

`-D with_util=` {yes | no}

>  Install utility tools. The default of this option is yes.
>
>  In case of cross-compilation, the frm tool is not installed. On the other hand, in case of native compile environment, `-D with_util=yes` indicates that all utilities are installed at the same time as CDMlib is compiled.
>
>  If `with_NetCDF=no` is specified, NetCDF function will not be supported.
>
>  If `with_MPI=no`, fconv and upacs2dfi will not be built.
>  

`-D with_MPI=` {yes | no}

>  If you use an MPI library, specify `with_MPI=yes`, the default is yes.

`-D with_TP =` *TextParser_directory*

> Specify the directory path that TextParser is installed.

`-D with_CPM=` {*CPMlib_directory* | no}

> If you want to build supplied utility tools, specify the directory path that CPM library is installed. CDMlib itself doesn't need CPMlib, but fconv depends on. If not specified this option, fconv is not compiled.

`-D with_NetCDF=` {no | *NetCDF_directory*}

> Specify the directory path that netCDF library is installed when you use netCDF4 file format. See the section of NetCDF support.

`-D with_HDF=installed_directory`
> Specify the directory path that HDF library is installed when you use HDF file format. See the section of HDF5 support.

`-D with_example=` {no | yes}

>  This option turns on compiling sample codes. The default is no.

`-D enable_BUFFER_SIZE=` {no | *size*}

>  Specify read/write buffer size in MB for accelerating file I/O performance. The default is no.

The default compiler options are described in `cmake/CompilerOptionSelector.cmake` file. See BUILD OPTION section in CMakeLists.txt in detail.


## Configure Examples

`$ export CDM_HOME=hogehoge`

In following exsmples, assuming that TextParser, CPMlib, and NetCDF library are installed under the CDM_HOME directory. If not, please specify applicable directory paths.

### INTEL/GNU compiler

~~~
$ cmake -DINSTALL_DIR=${CDM_HOME}/CDMlib -Dwith_MPI=yes -Dwith_util=yes -Dwith_example=yes -Dwith_TP=${CDM_HOME}/TextParser -Dwith_CPM=${CDM_HOME}/CPMlib -Dwith_HDF=no -Dwith_NetCDF=no -Denable_BUFFER_SIZE=no ..
~~~


### FUJITSU compiler / FX10, FX100, K on login nodes (Cross compilation) and Fujitsu TCS environment for intel PC

~~~
$ cmake -DINSTALL_DIR=${CDM_HOME}/CDMlib \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx10.cmake \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dwith_util=yes \
            -Dwith_TP=${CDM_HOME}/TextParser \
            -Dwith_CPM=${CDM_HOME}/CPMlib \
            -Dwith_HDF=no \
            -Dwith_NetCDF=no \
            -Denable_BUFFER_SIZE=no ..

$ cmake -DINSTALL_DIR=${CDM_HOME}/CDMlib \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx100.cmake \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dwith_util=yes \
            -Dwith_TP=${CDM_HOME}/TextParser \
            -Dwith_CPM=${CDM_HOME}/CPMlib \
            -Dwith_HDF=no \
            -Dwith_NetCDF=no \
            -Denable_BUFFER_SIZE=no ..

$ cmake -DINSTALL_DIR=${CDM_HOME}/CDMlib \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dwith_util=yes \
            -Dwith_TP=${CDM_HOME}/TextParser \
            -Dwith_CPM=${CDM_HOME}/CPMlib \
            -Dwith_HDF=no \
            -Dwith_NetCDF=no \
            -Denable_BUFFER_SIZE=no ..

$ cmake -DINSTALL_DIR=${CDM_HOME}/CDMlib \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_intel_F_TCS.cmake \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dwith_util=yes \
            -Dwith_TP=${CDM_HOME}/TextParser \
            -Dwith_CPM=${CDM_HOME}/CPMlib \
            -Dwith_HDF=no \
            -Dwith_NetCDF=no \
            -Denable_BUFFER_SIZE=no ..
~~~


##### Note
- On Fujitsu machines(fx10, K, fx100), confirm appropriate directrory path for compiler environment.
- Before building, execute following command for clean. `$ make distclean`



### Compilation of staging tool(frm)

When you install staging tool onto a login node with cross-compiling environment, the tool must be compiled by a native GNU compiler on the login node. If the front-end login node does not have MPI library, specify `-Dwith_MPI=no` option.



### NetCDF support
Before compiling NetCDF, we must build HDF5. HDF5 allows to use the options of several compression, zip, szip, etc. Since CDMlib manage multiple files, NetCDF and HDF5 may compile in serial version.

SZIP:
~~~
$ cat config_szip.sh
#!/bin/sh
export CC=gcc
export CFLAGS=-O3
#
./configure --prefix=$1

$ cd szip-x.x
$ config_szip.sh ${FFV_HOME}/SZIP
$ make
$ sudo make install
~~~

HDF5:
~~~
$ cat config_hdf5.sh
#!/bin/sh
export CC=gcc
export CFLAGS=-O3
export FC=gfortran
export FCFLAGS=-O3
#
./configure --prefix=$1 --enable-fortran --enable-cxx --with-zlib=/usr/bin --with-szip=${FFV_HOME}/SZIP

$ cd hdf5-x.x.x
$ config_hdf5.sh ${FFV_HOME}/HDF5
$ make
$ sudo make install
~~~

NetCDF:
~~~
$ cat config_netcdf.sh
#!/bin/sh
export CC=gcc
export CFLAGS=-O3
export CXX=g++
export CXXFLAGS=-O3
export FC=gfortran
export FCFLAGS=-O3
export CPPFLAGS=-I${FFV_HOME}/include
export LDFLAGS=-L${FFV_HOME}/lib
#
./configure --prefix=$1

$ cd netcdf-x.x.x
$ config_netcdf.sh ${FFV_HOME}/NetCDF
$ make
$ sudo make install
~~~


## EXAMPLES

* If you specify the test option by `-Denable_example=yes`, you can
execute the intrinsic tests by;

	`$ make test` or `$ ctest`

* The detailed results are written in `BUILD/Testing/Temporary/LastTest.log` file.
Meanwhile, the summary is displayed for stdout.




## CONTRIBUTORS

* Kenji     Ono      _keno@{cc.kyushu-u.ac, riken}.jp_
* Yasuhiro  Kawashima
* Eri       Takebe
* Yoichi    Tanaka
* Syoyo     Fujita
* Jorji     Nonaka
