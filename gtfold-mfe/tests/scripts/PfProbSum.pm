#!/usr/bin/perl
package PfProbSum;
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
  my $rnadir = $Config{"G_RNAFOLD_DIR"};
  my $workdir = $Config{"G_WORK_DIR"};
  my $paramdir = $gtdir."/../../rna-scoring/";

  my $key;
  my $value;
  my %new_hash = (%local_sequences, %Sequences);

  my $b2ct = "$workdir/b2ct";
  my $b2ct_compile = "gcc -o $b2ct ../cprogs/b2ct.c";
  my $b2ct_result = system("$b2ct_compile");

  if ($b2ct_result ne 0) {
    $logger->error("TEST_FAILED: Could not compile b2ct\n");
	return;
  }

  while (($key, $value) = each(%new_hash)) {

	  my $seqname = $key;
	  my $seqfile = $value;

  	  my $U1 = getPFvalue($seqfile, $gtdir);
	  #print($U1."\n");

  	  my $rnacmd;
	  my $subopt_output;
	  $rnacmd  = "$rnadir/RNAsubopt -s DAT -e 20 < $seqfile | grep \"[().]\"";
	  $subopt_output = `$rnacmd`;

	  $subopt_output =~ s/-.*//g;

      my @structures = split(/\n/,$subopt_output);
	  #print (scalar(@structures)."\n");

	  my $sum = 0;
	  my $structure;
	  foreach $structure (@structures) {
   	   	 my $b2ct_input = `cat $seqfile | grep "[ACGU]"`."$structure"."(1)";
	  	 #print "B2CT Input File:\n $b2ct_input\n";

	     my $ctFilePath = $workdir."/$seqname.ct";
		 my $b2ct_cmd = "echo \"$b2ct_input\" | $b2ct > $ctFilePath";
		 #print("$ctFilePath\n");
		 system("$b2ct_cmd"); 
		 #print ("here\n");
		 #print ("$ctFilePath");
		 #system("cat $ctFilePath");
         my $energy = getDSscore($ctFilePath, $paramdir);
		 #print "here1 $energy\n";
		 #my $x = $energy + 1;
		 #print $x."\n";
		 chomp($energy);
		 #exit(1);
		 #print $U1."\n";
         my $prob = getProbability($energy, $U1);
         $sum = $sum + $prob;
		 #print $prob."\n";
	  }

	  print("Final Sum: $sum\n");
	  if ($sum >= 0.95 || $sum <= 1.05) {
		$logger->info("TEST PASSED: Probability sum of all structures for $seqname is $sum");
	  }
	  else {
		$logger->info("TEST FAILED: Probability sum of all structures for $seqname is $sum");
	  }
  }
}

sub getProbability{
 my($e, $U)=@_;
 #P(I) = exp[-E(S, I)/RT]/U1;
 my $val = exp(-$e / 0.616333181);
 #print $val."\n";
 return $val/$U;
}

sub getPFvalue{
 my($seq, $gtdir)=@_;
#./gtboltzmann --partition pfTestSeqDB/combseq1/combseq1.seq
 my $cmd = "$gtdir/gtboltzmann --partition -dS".$seq;
 my $output = `$cmd`;
# print("output is:\n\n\n".$output);
my @lines = split(/\n/, $output);
#print($lines[(scalar@lines)-1]);
#print("\n\n");
return $lines[(scalar@lines)-2];
}

sub getDSscore{
 my($ctFilePath1, $paramDir1)=@_;
#./RNAScoring --dS ct_file_path 
 my $cmd = "$paramDir1/RNAScoring --dS --param-dir ".$paramDir1." ".$ctFilePath1;
 #print($cmd);
 my $output = `$cmd`;
	#printf("output:\n$output\n");
    my $lastline = `echo "$output" | grep \"Tree score is\" | sed 's/[A-Za-z ]*//'`;
    #print "here2 $lastline \n";
return $lastline;
}
1;
