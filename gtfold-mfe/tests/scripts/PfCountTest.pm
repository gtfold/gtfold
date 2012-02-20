#!/usr/bin/perl
package PfCountTest;
use strict;
use warnings;
use File::Basename;

sub test()
{
  my(%Config) = %{$_[1]};
  my(%Sequences) = %{$_[2]};
  my(%local_sequences) = %{$_[3]};
  my $logger = $_[4];

  my $gtdir = $Config{"G_GTFOLD_DIR"};
  my $unadir = $Config{"G_UNAFOLD_DIR"};
  my $workdir = $Config{"G_WORK_DIR"};

  my $key;
  my $value;
  my %new_hash = %local_sequences;

  while (($key, $value) = each(%new_hash)) {

	  my $seqname=$key;
	  my $seqfile = $value;
	  my $gtcmd;
      my $dirname = dirname($seqfile);
   	  my $count_file = "$dirname/$seqname.pfcount";
	  $gtcmd  = "$gtdir/gtboltzmann --pfcount $seqfile";
	  my $gtcount = `$gtcmd | tail -n 1`;
	  my $expected_count = `cat $count_file`;
	  if ($gtcount eq $expected_count) {
		$logger->info("TEST PASSED: $seqname: Counts matched : $gtcount");	
	  }
	  else {
		$logger->info("TEST FAILED: $seqname: Counts did not match. GTFOLD returned $gtcount, expected $expected_count.");
	  }

  }
}
1;
