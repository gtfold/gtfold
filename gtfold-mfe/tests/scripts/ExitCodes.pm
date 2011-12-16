#!/usr/bin/perl
package ExitCodes;
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
  my %new_hash = (%local_sequences, %Sequences);
  # For testing pseudoknot detection logic
  # we may only test the sequences specified in local_sequences
  while (($key, $value) = each(%new_hash)) {

	  my $seqname=$key;
	  my $path;
	  my $suffix;

	  my $seqfile = $value;
    my $dirname = dirname($seqfile);
    my $constraint_file = "$dirname/$seqname.constraint";
	  my $gtout  = "$workdir$seqname-gt";
	  my $gtoutfilename  = $workdir."$seqname-gt.ct";

	  my $gtcmd;
    if (-e $constraint_file) {
      $gtcmd  = "$gtdir/gtfold -c $constraint_file -m $seqfile -o $gtout > /dev/null 2>&1";
    }
    else {
      $gtcmd  = "$gtdir/gtfold -m $seqfile -o $gtout > /dev/null 2>&1";
    }

	  my $x = system("$gtcmd") >> 8;

    my $expected_result_file = "$dirname/$seqname.expectedresult";
    my $expected_result;
    if (-e $expected_result_file) {
      $expected_result = `head -n 1 $expected_result_file`;     
    }
    else {
      $expected_result = 0;
    }

    if ($expected_result =~ $x) {
      $logger->info("TEST PASSED: $seqname: Return Value $x matched expected output $expected_result");
    }
    else {
      $logger->error("TEST FAILED: $seqname: Return Value $x did not match expected output $expected_result");
    }
  }
}
1;
