#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use autodie;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use List::Util qw(sum);
use Data::Dumper;

################################################################################
# Option parsing
################################################################################
my %opts;
GetOptions( \%opts, "help|h", "man|m", "reference|ref|r:s", "calib|c:s",
    "excel|x", "cterr|e", "sep|s=s" )
    or pod2usage(2);

pod2usage(1) if $opts{"help"};
pod2usage( -exitstatus => 0, -verbose => 2 ) if $opts{"man"};

################################################################################
# Prototyping
################################################################################
sub calc_stats(@);
sub filter_outliers($$@);

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

    my $gene = $fields[ $sds_cols{detector} ];

    my $sample = $fields[ $sds_cols{sample} ];

    my $ct = $fields[ $sds_cols{ct} ];

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
    for my $sample ( keys %{ $qpcr_readings{$gene} } ) {
        my @ct_vals = @{ $qpcr_readings{$gene}->{$sample}->{ct} };

        my ( $mean, $stddev ) = calc_stats(@ct_vals);

        #warn "Old ct_vals: @ct_vals\n";

        # TODO When implemented, the --cterr option would
        # replace $stddev, if it is supplied. Otherwise, old
        # $stddev is used.
        warn "Filtering $gene, sample $sample for outlying values\n";
        ( $mean, $stddev, @ct_vals )
            = filter_outliers( $mean, $stddev, @ct_vals );
        warn "\n\n";

        #warn "New ct_vals: @ct_vals\n";

        $qpcr_readings{$gene}->{$sample}->{ct} = \@ct_vals;

        $qpcr_readings{$gene}->{$sample}->{mean} = $mean;
    }
}

# Now, all Ct values have been scanned for significant outliers, and the mean
# has been calculated. Time to get the delta_Ct, delta_delta_Ct and relative expression
my @genes = keys %qpcr_readings;

# Print a header for the output
my $sep = $opts{sep} // "\t";
say join( $sep,
    qw(Gene_name Sample_name Ct_values Mean_Ct delta_Ct delta_delta_Ct Relative_expression Comment)
);

# Start at "first gene" in list
for my $gene ( 0 .. @genes - 1 ) {

    # Go through all its samples and use the current sample as the reference
    my $gene_name = $genes[$gene];
    my @samples   = keys %{ $qpcr_readings{$gene_name} };

    for my $ref_sample ( 0 .. @samples - 1 ) {
        my $sample_name = $samples[$ref_sample];

      # Since this is the reference, fetch its mean and set its delta_Ct,
      # dd_Ct and re_express to 0, 0 and 1, respectively
      # TODO This needs to change for when manually setting reference
      # FIXME VERY cluttered because of nested structures
        my $ref_Ct   = $qpcr_readings{$gene_name}->{$sample_name}->{ct};
        my $ref_mean = $qpcr_readings{$gene_name}->{$sample_name}->{mean};

        my ( $ref_delta_Ct, $ref_d_delta_Ct, $ref_rel_express ) = ( 0, 0, 1 );

        say join(
            $sep,
            (   $gene_name,       $sample_name,
                "NEEDIMPL",       $ref_mean,
                $ref_delta_Ct,    $ref_d_delta_Ct,
                $ref_rel_express, "When used as reference"
            )
        );

        #Now use its mean for each subsequent sample
        for my $target_sample ( 0 .. @samples - 1 ) {
            next if ( $target_sample == $ref_sample );
            my $t_sample_name = $samples[$target_sample];
            my $t_Ct
                = $qpcr_readings{$gene_name}->{ $samples[$target_sample] }
                ->{ct};
            my $t_mean
                = $qpcr_readings{$gene_name}->{ $samples[$target_sample] }
                ->{mean};
            my $t_delta_Ct    = $t_mean - $ref_mean;
            my $t_d_delta_Ct  = $t_delta_Ct - $ref_delta_Ct;
            my $t_rel_express = 2**( -$t_d_delta_Ct );

            say join(
                $sep,
                (   $gene_name,     $t_sample_name,
                    "NEEDIMPL",     $t_mean,
                    $t_delta_Ct,    $t_d_delta_Ct,
                    $t_rel_express, "When compared to $sample_name"
                )
            );
        }
        say "";
    }

    # Then compare with other genes, if they exist.
}

#print Dumper(\%qpcr_readings);

################################################################################
# Subroutines
################################################################################
sub calc_stats(@) {
    my @ct_vals = @_;

    #TODO Probably better use Statistics::Basic
    my $sum  = sum(@ct_vals);
    my $mean = $sum / @ct_vals;

    my $diff_squares = 0;
    for my $val (@ct_vals) {
        $diff_squares += ( $val - $mean )**2;
    }

    # Use sample standard deviation (with Bessel's correction)
    # as R, Excel and others do
    my $stddev = sqrt( $diff_squares / ( @ct_vals - 1 ) );

    #warn "Sum: $sum, mean: $mean, stddev: $stddev\n";

    return ( $mean, $stddev );
}

sub filter_outliers($$@) {
    my ( $mean, $stddev, @values ) = @_;

    my @filtered;
    for my $val (@values) {
        if ( abs( $val - $mean ) <= $stddev ) {
            push @filtered, $val;
        }
        else {
            warn
                "Value $val is more than one standard deviation from sample mean: discarding...\n";
        }
    }

    ( $mean, $stddev ) = calc_stats(@filtered);

    return $mean, $stddev, @filtered;
}

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


