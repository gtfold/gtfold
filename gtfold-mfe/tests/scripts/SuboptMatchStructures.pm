#!/usr/bin/perl
package SuboptMatchStructures;
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
  my $rnadir = $Config{"G_RNAFOLD_DIR"};
  my $workdir = $Config{"G_WORK_DIR"};
  my $min_energy = $Config{"L_SUBOPTMATCHSTRUCTURES_MIN_ENERGY"};
  my $max_energy = $Config{"L_SUBOPTMATCHSTRUCTURES_MAX_ENERGY"};
  my $energy_stride = $Config{"L_SUBOPTMATCHSTRUCTURES_ENERGY_STRIDE"};

  if (not(defined($min_energy))) {
    $min_energy = 0;
  }
  if (not(defined($max_energy))) {
    $max_energy = 3;
  }
  if (not(defined($energy_stride))) {
    $energy_stride = 1;
  }

  my $key;
  my $value;
  my %new_hash = (%local_sequences);

  while (($key, $value) = each(%new_hash)) {

	  my $seqname=$key;
	  my $path;
	  my $suffix;

	  my $seqfile = $value;
    my $dirname = dirname($seqfile);
	  my $gtout  = "$workdir$seqname-gt";
    my $rnaout = "$workdir$seqname-rna";
	  my $gtoutfilename  = $workdir."$seqname-gt.ct";

    my %gtstruct_hash;
    my %rnastruct_hash;

    my $energy;
    for ($energy=$min_energy; $energy<=$max_energy; $energy=$energy+$energy_stride) {

      my $gtfile = $gtout.$energy;
      my $rnafile = $rnaout.$energy."_ss.txt";
  	  my $gtcmd;
      my $rnacmd;
      $gtcmd  = "$gtdir/gtsubopt --subopt $energy -o $gtfile $seqfile > /dev/null 2>&1";
      $rnacmd  = "$rnadir/RNAsubopt -s DAT -e $energy < $seqfile > $rnafile";

	    system("$gtcmd");
	    system("$rnacmd");

      my $gtsubopt_file = $gtfile."_ss.txt";
      my $rnasubopt_file = $rnafile;

      my @gtstructs = split( '[ \n][ \n]*', `cat $gtsubopt_file | sed '/[A-Z]/d' | sed '/*\$/d'`);
      my @rnastructs = split( '[ \n][ \n]*', `cat $rnasubopt_file | sed '/[A-Z]/d' | sed '/*\$/d'`);

	  my $num_elements = scalar(@gtstructs);
      if ( $num_elements % 2  == 1) {
		@gtstructs = @gtstructs[1,$num_elements-1];
	  }
   	  $num_elements = scalar(@rnastructs);
      if ( $num_elements % 2  == 1) {
		@rnastructs = @rnastructs[1,$num_elements-1];
	  }
      %gtstruct_hash = @gtstructs;
      %rnastruct_hash = @rnastructs;

      my $diff_file_name = "$workdir$seqname\_$energy.diff"; 

      open (DIFF_FILE, ">> $diff_file_name");

      my $diff = 0;
      foreach my $key ( keys %gtstruct_hash )
      {
        my $change = $gtstruct_hash{$key} 
        unless ( (exists $rnastruct_hash{$key}) && $gtstruct_hash{$key} eq $rnastruct_hash{$key} );
        # gotta check for existence to quiet warnings.
        if(defined($change) && $change ne "") {
          print DIFF_FILE "Key is ".$key. "\nGTFOLD Energy: ". $gtstruct_hash{$key};
          if (exists $rnastruct_hash{$key}) {
            print DIFF_FILE "\nRNASubopt Energy: ". $rnastruct_hash{$key};
          }
          else {
            print DIFF_FILE "\nRNASubopt Energy: No corresponding structure";
          }
          print DIFF_FILE "\n";
          $diff = 1;
        }
      }

      foreach my $key ( keys %rnastruct_hash )
      {
        my $change = $rnastruct_hash{$key} 
        unless ( (exists $gtstruct_hash{$key}) && $gtstruct_hash{$key} eq $rnastruct_hash{$key} );
        # gotta check for existence to quiet warnings.
        if(defined($change) && $change ne "") {
          print DIFF_FILE $key. "\nRNASubopt Energy: ". $rnastruct_hash{$key};
          if (exists $gtstruct_hash{$key}) {
            print DIFF_FILE "\nGTFOLD Energy: ". $gtstruct_hash{$key};
          }
          else {
            print DIFF_FILE "\nGTFOLD Energy: No corresponding structure";
          }
          print DIFF_FILE "\n";
          $diff = 1;
        }
      }

      close DIFF_FILE;

      if ($diff == 0) {
        $logger->info("TEST PASSED: $seqname: energy delta = $energy: Suboptimal Structures matched");
      }
      else {
        $logger->error("TEST FAILED: $seqname: energy delta = $energy: Suboptimal Structures not matched. Difference Written to file $diff_file_name");
      }
    }
  }
}
1;
