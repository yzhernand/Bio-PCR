
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

#use Bio::PCR::Well;
#use Bio::PCR::Sample;
#use Bio::PCR::Experiment;

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
   -file   : The name of the input file
   -format : Specify the format of the file.  Supported formats:

     sds22              SDS 2.2 format (default)

=cut

sub new {
    my ( $caller, @args ) = @_;
    my $self = {};

    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param;    # lowercase keys

    my $sdsfile = $param{'-file'};
    open my $fh, "<", $sdsfile;

    while ( my $line = <$fh> ) {
        next unless $line =~ /^\d/;

        #last if $line =~ /^Slope/;

        my @fields = split( /\t/, $line );

        my $gene = $fields[ $_sds_cols{detector} ];

        my $sample = $fields[ $_sds_cols{sample} ];

        my $ct = $fields[ $_sds_cols{ct} ];

        #    $found_ref = 1
        #        if ( $sample eq $opts{ref} );

        # qpcr_readings contains hashrefs based on sample
        # Each of those hashrefs use gene names as keys
        # and themselves point to another hashref.
        # That last nested hash will contain data for each sample
        push @{ $self->{qpcr_readings}->{$sample}->{$gene}->{ct} }, $ct;
    }

    close $fh;

    bless ($self, $caller);
    return $self;
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

    return $self->{qpcr_readings};
}
