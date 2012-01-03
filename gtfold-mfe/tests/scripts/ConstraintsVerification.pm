#!/usr/bin/perl
package ConstraintsVerification;
use strict;
use warnings;
use File::Basename;

require 'test_utils.pl';

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
	  my $path;
	  my $suffix;

	  my $seqfile = $value;
    my $dirname = dirname($seqfile);
    my $constraint_file = "$workdir/$seqname.constraint";
	  my $gtout  = "$workdir$seqname-gt";
	  my $gtoutfilename  = $workdir."$seqname-gt.ct";
	  my $gtcmd;
    $gtcmd  = "$gtdir/gtmfe $seqfile -o $gtout > /dev/null 2>&1";
	my $exitcode = system("$gtcmd") >> 8;
	if ($exitcode != 0) {
      $logger->error("TEST FAILED: $seqname : Assertion Failed, Exit Code incorrect: $exitcode.Could not execute MFE Calculation");
      next;	
	}
    my $len = get_seq_len_from_ctfile($gtoutfilename);
    my $CONSFILE;
    open($CONSFILE, '>', $constraint_file) or die "Couldn't open $_";
    generate_constraints(1, $len, $CONSFILE);
    close $CONSFILE;

    my $gtout_constraints = $gtout."_constraints";
    my $gtout_constraints_filename = "$gtout_constraints.ct";
    $gtcmd  = "$gtdir/gtmfe -c $constraint_file $seqfile -o $gtout_constraints > /dev/null 2>&1";
	$exitcode = system("$gtcmd") >> 8;

    if ($exitcode != 0) {
      $logger->error("TEST FAILED: $seqname : Assertion Failed, Exit Code incorrect: $exitcode. Constraint File: $constraint_file");
      next;
    }

    my @violated_constraints = verify_constraints($gtout_constraints_filename, $constraint_file);
    if (@violated_constraints) {
      my @vc = join(',', @violated_constraints);
#      print scalar(@violated_constraints);
      $logger->error("TEST FAILED: $seqname : Constraints Violated: Input Constraint File $constraint_file. Violated Constraints:\n @vc");
    }
    else {
      $logger->info("TEST PASSED: $seqname : Input Constraint File $constraint_file");
    }
  }
}

sub generate_constraints
{
  my $i = shift;
  my $j = shift;
  my $CONSFILE = shift;

  return if ($j-$i< 4); 
  
  my $r = rand();
  print $CONSFILE "F $i $j 1\n" if $r <0.4;
  print $CONSFILE "P $i $j 1\n" if $r >0.4 && $r<0.8;
  print $CONSFILE "P $i 0 1\n" if $r >0.8;

  my $range = $j-$i;
  my $ii = int(rand($range)) + $i;
  $range = $j-$ii;
  my $jj = int(rand($range)) + $ii;

  generate_constraints($i+1,$ii-1, $CONSFILE);
  generate_constraints($jj+1,$j-1, $CONSFILE);
}

sub get_seq_len_from_ctfile
{
  my $filename = shift;
  open(IN, "<$filename") || die ("Could not open $filename");
  my $header  = <IN>; # ignore first line
  close IN;
  my ($len, $dG) = split(' ', $header);
  return $len;
}

sub verify_constraints {

  my $ctfilename = shift;
  my $cffilename = shift;

  open CT, $ctfilename || die("Could not open $ctfilename");
  open CF, $cffilename || die("Could not open $cffilename");

  my %ct_hash = ();
  my @rna;
  my $count=1;
  my $header  = <CT>; # ignore first line
  while (<CT>)
  {
    my @ff = split(' ', $_);
    $ct_hash{$ff[0]} = $ff[4];
    $rna[$count] = $ff[1];
    ++$count;
  }

  my @violated_constraints = ();
  while(<CF>) 
  {
    chomp;
    my @ff = split(' ', $_);
    if ($ff[0] eq 'P' && $ff[2] == 0) { # case SS
      if ($ct_hash{$ff[1]} != 0) {
        print $_."\n";
        push(@violated_constraints, $_);
      }
    }
    elsif ($ff[0] eq 'P' && $ct_hash{$ff[1]} == $ff[2] && $ct_hash{$ff[2]} == $ff[1]) { # case P
        print $_."\n";
        print "$ct_hash{$ff[1]} $ct_hash{$ff[2]}\n";
      push(@violated_constraints, $_);
    }
    elsif ($ff[0] eq 'F' && $ff[2] == 0)	{
      my $b = canPair($rna[$ff[1]], $rna[$ff[2]]); 	
      if ($b == 1) {
        print $_."\n";
        push(@violated_constraints, $_);
      }
    }
  }

  close CT;
  close CF;

  return @violated_constraints;
}

sub canPair {
	 my ( $i, $j) = @_;
	# (A,U), (U,A), (C,G), (G,C), (G,U) and (U,G). 
	 return 1 if ($i eq 'A' && $j eq 'U') ;
	 return 1 if ($i eq 'U' && $j eq 'A') ;
	 return 1 if ($i eq 'C' && $j eq 'G') ;
	 return 1 if ($i eq 'G' && $j eq 'C') ;
	 return 1 if ($i eq 'G' && $j eq 'U') ;
	 return 1 if ($i eq 'U' && $j eq 'G') ;
	 return 0;
}
1;
