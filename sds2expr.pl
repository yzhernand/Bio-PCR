#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use GetOpt::Long qw(:config gnu_getopt);
use Pod::Usage;

################################################################################
# Option parsing
################################################################################
my %opts;
GetOptions(
        \%opts,
        "help|h",
        "man|m",
        "reference|ref|r:s",
        "calib|c:s",
        "excel|x"
          ) or pod2usage(2);

pod2usage(1) if $opts{"help"};
pod2usage( -exitstatus => 0, -verbose => 2 ) if $opts{"man"};

################################################################################
# Main
################################################################################

die "Usage: $0 <input SDS 2.2>"
    if @ARGV == 0;

my $filename = shift;

# Look up table for convenience. Has column indicies for needed data
my %col_positions = (Pos => 0, Sample => 1, Detector => 2, Ct => 5);

open my $fh, "<", $filename;

while (my $line = <$fh>) {
    next unless $line =~ /^\d/;
    last if $line =~ /^Slope/;

    my @fields = split(/\t/, $line);

    for my $col (@colnames) {
        if (!exists($new_col_positions{$col})) {
            print "\t";
            next;
        }

        print $fields[$new_col_positions{$col}], "\t";
    }
    print "\n";
}

close $fh;

################# POD Documentation ############

__END__

=head1 NAME

aln-manipulations.pl - Alignment tools based on BioPerl

=head1 SYNOPSIS

B<aln-manipulations.pl> [options] [alignment file]

=head1 DESCRIPTION

B<aln-manipulations.pl> will read an alignment file and 
do slice, display match characters, etc. By default, B<aln-manipulations.pl>
will assume both the input and the output are in clustalw format.

=head1 OPTIONS

=over 4

=item B<--help, -h>

Print a brief help message and exits.

=item B<--man, -m>

Print the manual page and exits.

=item B<--reference, --ref, -m>

Set the name of the reference/control sample.

This argument is optional. If it is not provided, sds2expr will simply go
through every sample and use every other sample as the reference.

=item B<--calib, -c>

Set the name of the calibrator sample.

This argument is optional. If it is not provided, the default is to use
the same sample name as --ref. If that is also not provided, then the
behavior is the same as described above, with the same sample used as
both reference and calibrator.

=item B<--excel, -x>

If specified, sds2expr will write an excel file to the same location of,
and with the same name as, the input file (with a .xls extension).

=back

=head1 EXAMPLES

Section under construction...

=head1 REQUIRES

  Perl 5.100

=head1 SEE ALSO

  perl(1)

=head1 AUTHOR

Yözen Hernández yzhernand at gmail dot com

=cut

##################### End ##########################


