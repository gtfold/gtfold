#!/usr/bin/perl
package SampleTest;
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

  my $key;
  my $value;
  my %new_hash = (%local_sequences, %Sequences);

  while (($key, $value) = each(%new_hash)) {

	  my $seqname=$key;
  }
}
1;
