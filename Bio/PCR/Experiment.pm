
=head1 NAME

Bio::PCR::Experiment - qPCR Experiment object

=head1 SYNOPSIS

    # DO NOT use this module directly. A Bio::PCR::Experiment object is created
    # by using Bio::PCR::PCRIO to open a file

    # Get a Bio::PCR::Experiment object somehow
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

A Bio::PCR::Experiment object is a collection of gene samples in a qPCR
trial. Samples are stored as Bio::PCR::Sample object, which in turn
contain information on qPCR wells in Bio::PCR::Well objects.

This module is probably better called 'Trial'.

=head1 FEEDBACK

=head1 AUTHOR - Yozen Hernandez

email - yzhernand at gmail dot com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::PCR::Well;
use strict;
use 5.010;

#use Bio::PCR::Sample;
#use Bio::PCR::Experiment;

# Used if a well does not have a position.
my $undef_sample_count = 0;

=head2 new

 Title   : new
 Usage   : my $experiment = Bio::PCR::Experiment->new(-name => $sample_name, -samples => \@sample_list);
 Function: Builds a new Bio::PCR::Experiment object
 Returns : Bio::PCR::Experiment
 Args    : a hash.  useful keys:
   -name    : Name of the sample or experiment
   -samples : Hash of gene names as keys and Bio::PCR::Sample objects as values
   -cterr   : Allowed error in Ct value for the experiment

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = {};

    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param;    # lowercase keys

    my $sample_name = $param{'-name'}
        // "Unnamed_experiment_" . $class::undef_sample_count++;

    my $samples = $param{'-samples'};
    my $ref     = $param{'-ref'} // undef;
    my $calib   = $param{'-calib'} // undef;

    $self->{SAMPLE_NAME} = $sample_name;
    $self->{SAMPLES}     = $samples;
    $self->{REF}         = $ref;
    $self->{CALIB}       = $calib;

    bless( $self, $class );
    return $self;
}

#=head2
#
# Title   :
# Usage   :
# Function:
# Returns :
# Args    :
#
#=cut
#
#sub {
#    my $self = shift;
#
#}
#
#=head2
#
# Title   :
# Usage   :
# Function:
# Returns :
# Args    :
#
#=cut
#
#sub {
#    my $self = shift;
#
#}

1;

