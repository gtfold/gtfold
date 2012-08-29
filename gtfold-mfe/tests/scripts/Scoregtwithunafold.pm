#!/usr/bin/perl
package Scoregtwithunafold;
use strict;
use warnings;
use File::Basename;

sub test()
{
  my(%Config) = %{$_[1]};
  my(%Sequences) = %{$_[2]};
  my $logger = $_[4];

  my $gtdir = $Config{"G_GTFOLD_DIR"};
  my $unadir = $Config{"G_UNAFOLD_DIR"};
  my $workdir = $Config{"G_WORK_DIR"};

  my $key;
  my $value;

  while (($key, $value) = each(%Sequences)) {

	  my $seqname=$key;
	  my $path;
	  my $suffix;

	  my $seqfile = $value;
	  my $gtout  = "$workdir$seqname-gt";
	  my $unaout = "$workdir$seqname-una";
	  my $gtoutfilename  = $workdir."$seqname-gt.ct";
	  my $unaoutfilename = $workdir."$seqname-una.ct";

	  my $gtcmd = "$gtdir/gtmfe -m $seqfile --unafold -o $gtout > /dev/null";
	  my $unacmd = "$unadir/hybrid-ss-min -s DAT $seqfile -o $unaout > /dev/null";
	
	  system("$gtcmd");
	  system("$unacmd");
	  my $gtfold_energy = `head -1 $gtoutfilename`;
	  my $unafold_energy = `head -1 $unaoutfilename`;

    if ($gtfold_energy =~ m/([0-9.-]+$)/) {
      $gtfold_energy = $1;
    }
    if ($unafold_energy =~ m/dG = ([0-9.-]+)/) {
      $unafold_energy = $1;
    }

    if ($gtfold_energy eq $unafold_energy) {
      $logger->info("TEST PASSED: $seqname: Energy Value matched for MFE -> $gtfold_energy");
    }
    else {
      $logger->error("TEST FAILED: $seqname: Energy Value Not Matched for MFE -> gtfold = $gtfold_energy, unafold = $unafold_energy");
    }
    #my $ctcmd1 = "/usr/local/bin/ct-energy $gtoutfilename";
	  #my $ctcmd2 = "/usr/local/bin/ct-energy -s DAT $unaoutfilename";
	  #system("$ctcmd1");
	  #print "ct-energy = ";
	  #system("$ctcmd2");
  }
}
1;
