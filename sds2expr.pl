#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use autodie;
use FindBin;                # Find the location of this script
use lib "$FindBin::Bin";    # to use it as a lib path
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use Bio::PCR::PCRIO::sds22;

################################################################################
# Option parsing
################################################################################
my %opts;
GetOptions(
    \%opts, "help|h", "man|m", "ref|reference|r:s",
    "calib|c:s" => \&opt_not_implemented,
    "excel|x"   => \&opt_not_implemented,
    "cterr|e=f", "sep|s=s", "use-gene|g"
) or pod2usage(2);

pod2usage(1) if $opts{"help"};
pod2usage( -exitstatus => 0, -verbose => 2 ) if $opts{"man"};

################################################################################
# Main
################################################################################

pod2usage(2) if ( @ARGV < 1 );

my $sep = $opts{'sep'};
my $ref = $opts{'ref'} // undef;
my $calib = $opts{'calib'} // undef;
my $cterr = $opts{'cterr'} // undef;
my $use_gene = $opts{'use-gene'} // undef;

my $filename = shift;

# Subject to change once PCRIO.pm is complete
my $sdsfile = Bio::PCR::PCRIO::sds22->new(-file => $filename, -sep => $sep, -ref => $ref, -cterr => $cterr, '-use-gene' => $use_gene);

my $experiments = $sdsfile->get_all_experiments;

for my $exp (@$experiments) {
    $exp->print_2ddct();
    print "\n";
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

=item B<--reference, --ref, -m> <name of sample>

Set the name of the reference/control sample.

This argument is optional. If it is not provided, sds2expr will simply go
through every sample and use every other sample as the reference.

=item B<--calib, -c> <name of sample>

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

=item B<--sep, -s> <separator>

Set the field delimiter in the output.

By default, all fields in output are tab-separated. By setting this option
to ',' (a comma) it is in effect the same as producing a CSV file (but on
standard output, not to a file).

Any string is a valid separator. The separator does not need to be only
one character.

=item B<--use-gene, -g>

If specified, will use the "Detector" column as the sample name, and the
"Sample" column as the "Detector" (the gene).

This is useful if you are running a qPCR on, for example, serial dilutions
of a single gene, where each dilution is your sample. By specifying this
option, the calculations will be made by comparing all dilutions of that
gene.

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


