#!/bin/bash


ffedir=$(pwd)
hdf5=hdf5-1.8.15-patch1
fftw=fftw-3.3.4


cat > Makefile.in <<EOF
# ----------------------------------------------------------
# Compiler configuration
# ----------------------------------------------------------
CC = cc
CCMM = \$(CC)
CFLAGS = -std=c99 -Wall -O3
CLIBS = -lm
AR = ar rcu
RANLIB = ranlib 


# ----------------------------------------------------------
# Optional dependencies (set to 1 if available)
# ----------------------------------------------------------
HAVE_MPI = 0
HAVE_HDF5 = 1
HAVE_FFTW = 1
HAVE_RNPL = 0


# ----------------------------------------------------------
# System library directories
# ----------------------------------------------------------
HDF5_HOME = $ffedir/$hdf5
FFTW_HOME = $ffedir/$fftw
RNPL_HOME = /usr/local
EOF


git clone https://github.com/jzrake/cow --branch=ffe-pspec
cd cow
ln -s ../Makefile.in .
cd ..


wget http://www.fftw.org/$fftw.tar.gz
tar -xvf $fftw.tar.gz
mv $fftw $fftw-src
cd $fftw-src
./configure --prefix=$ffedir/$fftw
make install
cd ..


wget http://www.hdfgroup.org/ftp/HDF5/current/src/$hdf5.tar.gz
tar -xvf $hdf5.tar.gz
mv $hdf5 $hdf5-src
cd $hdf5-src
./configure --prefix=$ffedir/$hdf5
make install
cd ..


make
