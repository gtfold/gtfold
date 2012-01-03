#!/usr/bin/perl
use strict;
use warnings;
use Module::Load;

use File::Basename;
use File::Path;
use Log::Log4perl qw(:easy);

require 'test_utils.pl';

my $configdir="../config";

Log::Log4perl::init( "$configdir/root-logger.conf" );
my $logger = Log::Log4perl->get_logger;


$logger->info("Starting Tests...");

# Read Parameter File
my $paramfile = "$configdir/test-params.conf";

# Create Work Directory
my $workdir = "../work";
$logger->info("Deleting directory $workdir");
rmtree($workdir, 0, 1);
$logger->info("Creating work directory $workdir");
mkdir $workdir;

my %Config;

%Config = load_config_file($paramfile);

###### Prepare a Hashmap of Sequences ######

my %Sequences;
my @seqdir_arr = @{$Config{"G_SEQUENCE_DIR"}};

my $seq_include_regex = $Config{"G_INCLUDE_SEQUENCES"};
my $seq_exclude_regex = $Config{"G_EXCLUDE_SEQUENCES"};

%Sequences = find_sequences_from_dir(
                seq_include_regex => $seq_include_regex,
                seq_exclude_regex => $seq_exclude_regex,
                seq_dirs => @seqdir_arr);

my $test_include_regex = $Config{"G_INCLUDE_TESTS"};
my $test_exclude_regex = $Config{"G_EXCLUDE_TESTS"};

my $test_list_file = $Config{"G_TEST_LIST_FILE"};

open(TESTLISTFILE, $test_list_file) || die("Could not open file: $test_list_file");

  while (<TESTLISTFILE>) {

      my $testname = $_;
      chomp($testname);

      $testname =~ s/^\s*//;     # Remove spaces at the start of the line
      $testname =~ s/\s*$//;     # Remove spaces at the end of the line
      if ( ($testname !~ /^#/) && ($testname ne "") ) {    # Ignore lines starting with # and blank lines

        my $test_include = (not defined($test_include_regex)) || ($testname =~ /$seq_include_regex/);
        my $test_exclude = (not defined($test_exclude_regex)) || ($testname !~ /$test_exclude_regex/);

        if ( $test_include && $test_exclude ) {
          my $uppertestname = uc($testname);
          my $local_sequence_dirname = "L_".$uppertestname."_SEQUENCE_DIR";
          my $local_sequence_include_regex = $Config{"L_".$uppertestname."_SEQUENCE_INCLUDE"};
          my $local_sequence_exclude_regex = $Config{"L_".$uppertestname."_SEQUENCE_EXCLUDE"};
          my @local_sequence_dir;
          my %local_sequences;

          if (defined($Config{$local_sequence_dirname})) {
            @local_sequence_dir = @{$Config{$local_sequence_dirname}};
            %local_sequences = find_sequences_from_dir(
                    seq_include_regex => $local_sequence_include_regex,
                    seq_exclude_regex => $local_sequence_exclude_regex,
                    seq_dirs => @local_sequence_dir);
          }

          my $module = $testname;
          load($module);
          $module->test(\%Config, \%Sequences, \%local_sequences, $logger);
        }
      }
   }


sub find_sequences_from_dir
{

#print @_[1]."\n";
#print @_[3];
#print $_[5i]."\n";
#my %args = @_;

my %seq_hash;

#if (defined $_[1]) {
	my $include_regex = $_[1];
#}
#if (defined $_[3]) {
	my $exclude_regex = $_[3];
#}
#if not(defined $_[5]) {
#	return %seq_hash;
#}

shift @_;
shift @_;
shift @_;
shift @_;
shift @_;

foreach (@_) {

  opendir(DIR, $_) || die $!;
  while (my $seqfile = readdir(DIR)) {

	  my $seqname;
	  my $path;
	  my $suffix;
	  ($seqname,$path,$suffix) = fileparse($seqfile, (".seq", ".constraint"));

	  if ((-d $seqfile) || $suffix ne '.seq') {
		  next;
	  }

    $seqfile = "$_$seqname$suffix";
    my $seq_include = (not defined($include_regex)) || ($seqname =~ /$include_regex/);
    my $seq_exclude = (not defined($exclude_regex)) || ($seqname !~ /$exclude_regex/);

    if ( $seq_include && $seq_exclude ) {
      $seq_hash{$seqname} = $seqfile;
      $logger->info("Selected Sequence ... $seqname");
    }
  }

}
  return %seq_hash;
}
