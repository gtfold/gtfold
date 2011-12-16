#!/bin/sh

BINDIR='.'
echo 'creating symlinks under' $BINDIR
ln -sf gtfold $BINDIR/gtmfe
ln -sf gtfold $BINDIR/gtsubopt
ln -sf gtfold $BINDIR/gtboltzmann

