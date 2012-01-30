#!/usr/bin/perl
use strict;
use warnings;

my $scoreCalcDir = "/home/users/msoni/Desktop/manoj_gatech/research/gtfold/scoring_code_shel/shelswenson-gtfold-cdc688c/rna-scoring";
my $seqDir = $scoreCalcDir."/d0d2dScomparison/sequences";
my $scoreCalcPrg = $scoreCalcDir."/RNAScoring";

my $easyoutfilepath = $seqDir."/easy/d0d2dScomparison_mfe_easy.txt";

my $easyseq1=$seqDir."/easy/biloba.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/easy/mobilis.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/easy/phytomonas.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/easy/spombe.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/easy/thermophila.ct";
generateResult($easyseq1,$easyoutfilepath);

sub generateResult{
 my($seq,$out)=@_;
 #my $str="Score is\n";
 #`echo $str`;
#print("$scoreCalcPrg $seq >> $out");
 `$scoreCalcPrg --d0 $seq >> $out`;
 `$scoreCalcPrg $seq >> $out`;
 `$scoreCalcPrg --dS $seq >> $out`;
 `$scoreCalcPrg --d2 $seq >> $out`;
}
