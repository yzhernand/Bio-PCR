#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use List::Util qw(sum);
use Data::Dumper;

################################################################################
# Option parsing
################################################################################
my %opts;
GetOptions( \%opts, "help|h", "man|m", "reference|ref|r:s", "calib|c:s",
    "excel|x", "cterr|e" )
  or pod2usage(2);

pod2usage(1) if $opts{"help"};
pod2usage( -exitstatus => 0, -verbose => 2 ) if $opts{"man"};

################################################################################
# Main
################################################################################

pod2usage(2) if ( @ARGV < 1 );

my $sdsfile = shift;

# Look up table for convenience. Has column indicies for needed data
my %sds_cols = ( pos => 0, sample => 1, detector => 2, ct => 5 );

# Hash of genes to a hashref of samples
my %qpcr_readings = ();

open my $fh, "<", $sdsfile;

while ( my $line = <$fh> ) {
    next unless $line =~ /^\d/;
    last if $line =~ /^Slope/;

    my @fields = split( /\t/, $line );

    my $gene   = $fields[ $sds_cols{detector} ];

    my $sample = $fields[ $sds_cols{sample} ];

    my $ct     = $fields[ $sds_cols{ct} ];

    # %qpcr_readings contains hashrefs based on gene
    # Each of those hashrefs use sample names as keys
    # and themselves point to another hashref.
    # That last nested hash will contain data for each sample
    push @{ $qpcr_readings{$gene}->{$sample}->{ct} }, $ct;
}

close $fh;

# By this point, all qPCR data will have been read.
#
# Begin filtering outliers from each sample, and
# compute the true average and stddev

for my $gene ( keys %qpcr_readings ) {
    for my $sample ( keys %{$qpcr_readings{$gene}} ) {
        my @ct_vals = @{ $qpcr_readings{$gene}->{$sample}->{ct} };

        my $sum = sum(@ct_vals);

        #TODO Probably better use Statistics::Basic
        my $mean = $sum/ @ct_vals;

        my $diff_squares = 0;
        for my $val ( @ct_vals ){
            $diff_squares += ($mean-$val)**2;
        }

        my $stddev = $diff_squares / @ct_vals;

        say "Sum: $sum, mean: $mean, stddev: $stddev";
    }
}

print Dumper(\%qpcr_readings);

################# POD Documentation ############

__END__

=head1 NAME

sds2expr.pl - Calculate relative expression change using the
2^(-Delta Delta Ct) method. Input must be an SDS 2.2.2 formatted text
file.

=head1 SYNOPSIS

B<sds2expr.pl> [options] <SDS 2.2.2 results text file>

=head1 DESCRIPTION

B<sds2expr.pl> takes as input a qPCR results file from SDS 2.2.2 and
calculates the relative expression change using the 2^(-Delta Delta Ct)
method.

The calculation can be done either using user-specified reference and
calibrator samples, or doing pairwise calculations were every sample
is used as a reference at least once.

Note that files produced by SDS 2.3 or higher are not supported. For
these, it is recommended to use the Bioconductor module HTqPCR for R.

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

=item B<--cterr, -e> <error (real number)>

Specify the amount of error allowed in Ct values. If at least one well
had a reading +/- "error" #TODO

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


