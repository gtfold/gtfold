#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

#die "Usage: ./gtmfe_test_options.pl <seq file>" if @ARGV < 1;
my $gtfoldBinDir = "/home/users/msoni/gtfold/src";
my $rnafoldBinDir = "/home/users/msoni/Desktop/manoj_gatech/research/viennaRNA/ViennaRNA-1.8.5/Progs";
my $unafoldBinDir = "/home/users/msoni/Desktop/manoj_gatech/research/unafold/unafold-3.8/scripts";

my @seqNames = ("/home/users/msoni/gtfold/src/pfTestSeqDB/combseq1/combseq1.seq", "/home/users/msoni/gtfold/src/pfTestSeqDB/combseq2/combseq2.seq", "/home/users/msoni/gtfold/src/pfTestSeqDB/combseq3/combseq3.seq", "/home/users/msoni/gtfold/src/pfTestSeqDB/combseq4/combseq4.seq", "/home/users/msoni/gtfold/src/pfTestSeqDB/seq600/seq600.seq", "/home/users/msoni/gtfold/src/pfTestSeqDB/ecoli16s_900/ecoli16s_900.seq");

sub runrnafoldtest(){
	my $seqName;
	foreach $seqName (@seqNames){
		rnafoldtest($seqName);
	}
}

sub rununafoldtest(){
	my $seqName;
	foreach $seqName (@seqNames){
		unafoldtest($seqName);
	}
}

sub rnafoldtest(){
	my ($seqName) = @_;
	rungtfold($seqName, "--rnafold");
	runrnafold($seqName);
	print("----------------------------------------------\n\n\n");
}

sub unafoldtest(){
	my ($seqName) = @_;
	rungtfold($seqName, "--unafold");
	rununafold($seqName);
	print("----------------------------------------------\n\n\n");
}



sub runrnafold(){
	my ($seqName) = @_;
	my $rnafoldCmd = "cat ".$seqName." | ".$rnafoldBinDir."/RNAfold";
	print("Executing $rnafoldCmd\n");
	print("----------------------------------------------\n");
	system($rnafoldCmd);
}

sub rununafold(){
	my ($seqName) = @_;
	my $cmd = $unafoldBinDir."/UNAFold.pl ".$seqName;
	print("Executing $cmd\n");
	print("----------------------------------------------\n");
	system($cmd);
}


sub rungtfold(){
	my ($seqName, $options) = @_;
	my $gtfoldCmd = $gtfoldBinDir."/gtmfe"." --silent ".$options." ".$seqName;
	print("Executing $gtfoldCmd\n");
	print("----------------------------------------------\n");
	system($gtfoldCmd);
}

#runrnafoldtest();
rununafoldtest();
print("gtmfe_test_options program completed successfully\n.")

