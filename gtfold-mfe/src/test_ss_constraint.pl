#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

die "Usage: ./test_ss_constraint.pl <seq file>" if @ARGV < 1;

my $seqName = $ARGV[0];

system("gcc -o test_ss_constraint test_ss_constraint.c");
system("./test_ss_constraint $seqName");
