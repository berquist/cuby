#!/bin/bash

LIBDIR=$HOME/cuby4_libs
mkdir $LIBDIR
cd $LIBDIR

# Download and unpack LAPACK
# Older version is used so that autoconf is available
wget http://www.netlib.org/lapack/lapack-3.4.0.tgz
tar -xf lapack-3.4.0.tgz
cd lapack-3.4.0

# Download autoconf and run it
wget http://users.wfu.edu/cottrell/lapack/lapack-3.4.0-autoconf.tar.gz
tar -xf lapack-3.4.0-autoconf.tar.gz
./configure --prefix=$LIBDIR --enable-single=no --enable-complex=no --enable-complex16=no --enable-shared

# Make and install
make
make install

# Message

echo "================================================================================"
echo "To use the libraries, add following to ~/.cuby4":
echo "extension_lib_dirs:"
echo "  - $LIBDIR"
echo "================================================================================"

