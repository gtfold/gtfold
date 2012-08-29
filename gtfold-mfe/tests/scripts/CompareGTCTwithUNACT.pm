#!/usr/bin/perl
package CompareGTCTwithUNACT;
use strict;
use warnings;
use Error;
use File::Basename;
use Exception::Class;

require "compare_ct.pl"; 

sub test()
{
  my(%Config) = %{$_[1]};
  my(%Sequences) = %{$_[2]};
  my $logger = $_[4];

  my $epsilon = 0.001;
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

    my @result;
    eval {
      @result = compare_ct_structures($gtoutfilename, $unaoutfilename);
    };
    my $exception;
    if ($exception = Exception::Class->caught('CompareCTException')) {
      $logger->error("TEST FAILED: $seqname: $exception\n");
      next;
    }

    my $diff1 = abs(1 - $result[0]);
    my $diff2 = abs(1 - $result[1]);
    if ($diff1 < $epsilon || $diff2 < $epsilon) {
      $logger->info("TEST PASSED: $seqname: Specificity = $result[0], Selectivity = $result[1]");
    }
    else {
      $logger->info("TEST FAILED: $seqname: Specificity = $result[0], Selectivity = $result[1]\n".
                    "More Information: $result[2] $result[3]");
    }
  }
}
1;
