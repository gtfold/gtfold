#!/usr/bin/perl
use strict;
use warnings;

my $scoreCalcDir = "/home/users/msoni/Desktop/manoj_gatech/research/gtfold/scoring_code_shel/shelswenson-gtfold-cdc688c/rna-scoring";
my $seqDir = $scoreCalcDir."/d0d2dScomparison/sequences";
my $scoreCalcPrg = $scoreCalcDir."/RNAScoring";

my $easyoutfilepath = $seqDir."/medium/d0d2dScomparison_mfe_medium.txt";

my $easyseq1=$seqDir."/easy/biloba.ct";
#--------------------------------------------
$easyseq1=$seqDir."/medium/japonica.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/medium/norvegicus.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/medium/sanguinea.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/medium/taurus.ct";
generateResult($easyseq1,$easyoutfilepath);

$easyseq1=$seqDir."/medium/waltl.ct";
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
