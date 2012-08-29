#!/bin/sh
if [ $# != 1 ]
then
echo "Usage: sh createsymlinks.sh <INSTALLDIR>"
exit -1
fi

BINDIR=$1 
echo 'creating symlinks under' $BINDIR
cd $BINDIR
ln -sf gtfold gtmfe
ln -sf gtfold gtsubopt
ln -sf gtfold gtboltzmann

