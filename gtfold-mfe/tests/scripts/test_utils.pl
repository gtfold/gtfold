#!/usr/bin/perl
use strict;
use File::Path;

sub load_config_file()
{
  my $paramfile = $_[0];
  my %Config;
  my $Name;
  my $Value;
  my $line;
  my $orig_value;
  open(PARAMFILE, $paramfile) || die("Could not open file: $paramfile");

  while (<PARAMFILE>) {

    my $line = $_;
    chomp($line);

    $line =~ s/^\s*//;     # Remove spaces at the start of the line
    $line =~ s/\s*$//;     # Remove spaces at the end of the line
    if ( ($line !~ /^#/) && ($line ne "") ) {    # Ignore lines starting with # and blank lines
      ($Name, $Value) = split (/=/, $line);    # Split each line into name value pairs

      if ($Name =~ /.*_INCLUDE_SEQUENCES|.*_EXCLUDE_SEQUENCES|G_INCLUDE_TESTS|G_EXCLUDE_TESTS/) {
        $orig_value= %Config->{$Name};
        $Value =~ s/[.]/\\./g;
        $Value =~ s/[*]/.\*/g;
        if ($orig_value ne "") {
          $Value = $orig_value."|".$Value;
        }
      }

      if ($Name =~ /.*_SEQUENCE_DIR/) {
        my @seqdir_arr = %Config->{$Name};
        push (@{%Config->{$Name}}, $Value);
        next;
      }

      %Config->{$Name} = $Value;     # Create a hash of the name value pairs

    }
  }

  $Value = %Config->{"G_INCLUDE_SEQUENCES"};

#  printf($Value."\n");
#  if ("rad" =~ /\G$Value/) {
#    printf("matched1\n");
#  }
#  printf($Value."\n");
#  if ("a.a.dvdsf" =~ /\G(rad|a\..*)/g) {
#    printf("matched2\n");
#  }
#  printf($Value."\n");
#  if ("daf" =~ /\G$Value/) {
#    printf("matched3\n");
#  }
#  printf($Value."\n");
#  if ("d.afdsf" =~ /\G$Value/) {
#    printf("matched4\n");
#  }
#  printf($Value."\n");
#  if ("acc1" =~ /\G$Value/) {
#    printf("matched5\n");
#  }
#  printf($Value."\n");
#  if ("acc11" =~ /\G$Value/) {
#    printf("matched6\n");
#  }
#  printf($Value."\n");
#  if ("a..afjlhfj" =~ /\G$Value/) {
#    printf("matched7\n");
#  }
#  printf($Value."\n");
#  if ("afdjs" =~ /\G$Value/) {
#    printf("matched8\n");
#  }

#  printf ("@{%Config->{'G_SEQUENCE_DIR'}}\n");

  return %Config;
}
1;
