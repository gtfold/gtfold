#!/usr/bin/perl
package StochasticEnergyTest;
use strict;
use warnings;

sub test()
{
  my(%Config) = %{$_[1]};
  my(%Sequences) = %{$_[2]};
  my(%local_sequences) = %{$_[3]};
  my $logger = $_[4];

  my $gtdir = $Config{"G_GTFOLD_DIR"};
  my $unadir = $Config{"G_UNAFOLD_DIR"};
  my $workdir = $Config{"G_WORK_DIR"};
  my $rnascoring_dir = "$gtdir/../../rna-scoring";

  my $key;
  my $value;
  my %new_hash = (%local_sequences, %Sequences);

  while (($key, $value) = each(%new_hash)) {

	  my $seqname=$key;
      my $seqfile = $value;

      my $summary_filename = "$seqname.summary";
      my $summaryfile = "$workdir/$summary_filename";
      my $gtcmd = "$gtdir/gtboltzmann --sample 2 --dump --dump_dir $workdir --dump_summary $summary_filename $seqfile";
      print $gtcmd."\n";
      system("$gtcmd");

      my $output_file = "$workdir/$seqname.stochastic.out";
      my $error_file = "$workdir/$seqname.stochastic.err";
      my $rnascoring_cmd = "$rnascoring_dir/RNAScoring --pf_test $summaryfile $output_file $error_file";

      system("$rnascoring_cmd");

	  if (-z "$error_file") {
		$logger->info("TEST PASSED: sequence $seqname");
	  }
	  else {
		$logger->info("TEST FAILED: sequence $seqname. For more information refer $output_file and $error_file");
	  }

      system("$gtcmd");
  }
}
1;
