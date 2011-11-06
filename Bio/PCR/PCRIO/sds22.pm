
=head1 NAME

Bio::PCR::PCRIO::sds22 - Driver for reading qPCR result files
produced by SDS 2.2.2

=head1 SYNOPSIS

    use strict;
    use 5.010;
    use Bio::PCR::PCRIO::sds22; # Will change when PCRIO.pm is ready
    
    my $filename = "qpcr_results.txt";
    
    my $sdsfile = Bio::PCR::PCRIO::sds22->new(-file => $filename);
    
    # Hash reference
    my $experiments = $sdsfile->get_all_experiments();
    
    for my $exp ( keys( %$experiments ) ){
        say $experiments->{$exp}; # Prints out all experimental data
    }

=head1 DESCRIPTION

This is (will be) the driver for reading SDS 2.2 formatted
files for PCR data

=head1 FEEDBACK

=head1 AUTHOR - Yozen Hernandez

email - yzhernand at gmail dot com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::PCR::PCRIO::sds22;
use strict;
use 5.010;
use Carp;

use Bio::PCR::Well;
use Bio::PCR::Sample;
use Bio::PCR::Experiment;
#use Data::Dumper;

our $VERSION = '0.5';

# Look up table for convenience. Has column indicies for needed data
my %_sds_cols = ( pos => 0, sample => 1, detector => 2, ct => 5 );

# Hash of genes to a hashref of samples
#my %qpcr_readings = ();

#my $found_ref     = 0;

=head2 new

 Title   : new
 Usage   : my $sdsfile = Bio::PCR::PCRIO::sds22->new(-file => filename);
 Function: Builds a new Bio::PCR::PCRIO object 
 Returns : Bio::PCR::PCRIO
 Args    : a hash.  useful keys:
   -file    : The name of the input file
   -cterr   : Allowed error in Ct value
   -ref     : Reference sequence
   -calib   : Calibrator sequence
   -format  : Specify the format of the file.  Supported formats:

     sds22              SDS 2.2 format (default)

=cut

sub new {
    my ( $caller, @args ) = @_;
    my $self = {};

    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param;    # lowercase keys

    # Optional arguments
    $self->{SEP}      = $param{'-sep'};
    $self->{CTERR}    = $param{'-cterr'} // undef;
    $self->{REFSEQ}   = $param{'-ref'}   // undef;
    $self->{CALIBSEQ} = $param{'-calib'} // undef;

    ($_sds_cols{sample}, $_sds_cols{detector}) = ($_sds_cols{detector}, $_sds_cols{sample})
        if $param{'-use-gene'};

    my $sdsfile = $param{'-file'};
    open my $fh, "<", $sdsfile;

    while ( my $line = <$fh> ) {
        next unless $line =~ /^\d/;

        #last if $line =~ /^Slope/;

        my @fields = split( /\t/, $line );

        my $pos = $fields[ $_sds_cols{pos} ];

        my $gene = $fields[ $_sds_cols{detector} ];

        my $sample = $fields[ $_sds_cols{sample} ];

        my $ct = $fields[ $_sds_cols{ct} ];

        #$found_ref = 1
        #    if ( $sample eq $opts{ref} );

        unless ( defined( $self->{qpcr_readings}->{$sample}->{$gene} ) ) {
            my $sample_obj = Bio::PCR::Sample->new(
                -name  => $gene,
                -cterr => $self->{CTERR}
            );
            $self->{qpcr_readings}->{$sample}->{$gene} = $sample_obj;
        }

        if ( $ct =~ /Undetermined/ ) {
            carp "WARNING: Well $pos is Undetermined. Ignoring...\n";
            next;
        }

        my $well = Bio::PCR::Well->new( -ct => $ct, -pos => $pos );
        $self->{qpcr_readings}->{$sample}->{$gene}->add_well($well);
    }

    close $fh;

    bless( $self, $caller );
    $self->_make_experiments();

    return $self;
}

=head2 _make_experiments

 Title   : _make_experiments
 Usage   : $self->_make_experiments();
 Function: Creates experiment objects using the internal hash built from file
 Returns : None
 Args    : None

=cut

sub _make_experiments {
    my $self = shift;

    for my $sample ( keys %{ $self->{qpcr_readings} } ) {
        my $experiment = Bio::PCR::Experiment->new(
            -name    => $sample,
            -samples => $self->{qpcr_readings}->{$sample},
            -sep     => $self->{SEP},
            -ref     => $self->{REFSEQ},
            -calib   => $self->{CALIBSEQ}
        );
        push( @{ $self->{EXPERIMENTS} }, $experiment );
    }
}

=head2 get_all_experiments

 Title   : get_all_experiments
 Usage   : my $experiments = Bio::PCR::PCRIO::sds22->get_all_experiments();
 Function: Retrieves all experiments in a file
 Returns : Hash reference
 Args    : None

=cut

sub get_all_experiments {
    my $self = shift;

    return $self->{EXPERIMENTS};
}

1;

