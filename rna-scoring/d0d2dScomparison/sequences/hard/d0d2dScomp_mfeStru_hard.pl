#!/usr/bin/perl
use strict;
use warnings;

my $scoreCalcDir = "/home/users/msoni/Desktop/manoj_gatech/research/gtfold/scoring_code_shel/shelswenson-gtfold-cdc688c/rna-scoring";
my $seqDir = $scoreCalcDir."/d0d2dScomparison/sequences";
my $scoreCalcPrg = $scoreCalcDir."/RNAScoring";

my $easyoutfilepath = $seqDir."/hard/d0d2dScomparison_mfe_hard.txt";

my $easyseq1=$seqDir."/easy/biloba.ct";
#------------------------------------------
$easyseq1=$seqDir."/hard/arboreum.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/hard/canicula.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/hard/ecoli.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/hard/hsapiens.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/hard/volcanii.ct";
generateResult($easyseq1,$easyoutfilepath);

sub generateResult{
 my($seq,$out)=@_;
 #my $str="Score is\n";
 #`echo $str`;
 `$scoreCalcPrg --d0 $seq >> $out`;
 `$scoreCalcPrg $seq >> $out`;
 `$scoreCalcPrg --dS $seq >> $out`;
 `$scoreCalcPrg --d2 $seq >> $out`;
}
