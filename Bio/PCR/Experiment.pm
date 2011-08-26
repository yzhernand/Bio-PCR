
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

package Bio::PCR::Experiment;
use strict;
use 5.010;
use Carp;
use Data::Dumper;
use Bio::PCR::Sample;
use Bio::PCR::Well;


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
        // "Unnamed_experiment_";

    my $samples = $param{'-samples'};

    croak "Error when setting samples. '-samples' must be a hash ref\n"
        unless(ref($samples) eq 'HASH');

    my $ref     = $param{'-ref'} // undef;
    my $calib   = $param{'-calib'} // undef;
    croak "Error: reference sample \"$ref\" not found. ",
        "Please check name is exactly as it appears in input file. Exiting...\n"
        if(defined($ref) && (!exists $samples->{$ref} ));

    $self->{SAMPLE_NAME} = $sample_name;
    $self->{SAMPLES}     = $samples;
    $self->{REF}         = $ref;
    $self->{CALIB}       = $calib;

    bless( $self, $class );
    $self->_discard_empty_samples();

    return $self;
}


=head2 _calc_2ddct

 Title   : _calc_2ddct
 Usage   : $self->_calc_2ddct() # INTERNAL METHOD
 Function: Calculates the 2^-ddCt value for a sample, given a reference sample
 Returns : None
 Args    : Reference smple, Target sample. Both are Bio::PCR::Sample objects

=cut

sub _calc_2ddct{
    my ($self, $ref_sample, $target_sample) = @_;

    my $r_Ct = $self->{SAMPLES}->{$ref_sample}->get_avg_ct();
    my $t_Ct = $self->{SAMPLES}->{$target_sample}->get_avg_ct();
    my $t_dct = $t_Ct - $r_Ct;
    my $t_dd_ct = $t_dct - $self->{SAMPLES}->{$ref_sample}->get_dd_ct();
    my $t_rel_express = 2**(-$t_dd_ct);
    $self->{SAMPLES}->{$target_sample}->set_delta_ct($t_dct);
    $self->{SAMPLES}->{$target_sample}->set_dd_ct($t_dd_ct);
    $self->{SAMPLES}->{$target_sample}->set_2ddct($t_rel_express);
}


=head2 _calc_2ddct_all

 Title   : _calc_2ddct_all
 Usage   : $self->_calc_2ddct_all() # INTERNAL METHOD
 Function: Calculates the 2^-ddCt value for all samples in an experiment
 Returns : None
 Args    : None

=cut

sub _calc_2ddct_all{
    my $self = shift;
    
    for my $reference( keys %{ $self->{SAMPLES} } ) {
        if (defined $self->{REF}) {
            next unless ($reference eq $self->{REF});
        }

        $self->{SAMPLES}->{$reference}->set_delta_ct(0);
        $self->{SAMPLES}->{$reference}->set_dd_ct(0);
        $self->{SAMPLES}->{$reference}->set_2ddct(1);

        for my $target ( keys %{ $self->{SAMPLES} }) {
            next if ($target eq $reference);
            
            $self->_calc_2ddct($reference,$target);
        }

        say Dumper($self->{SAMPLES});
    }

}

=head2 _discard_empty_samples

 Title   : _discard_empty_samples
 Usage   : $self->_discard_empty_samples # INTERNAL METHOD
 Function: Goes through all samples and discards any which have no methods
 Returns :
 Args    :

=cut

sub _discard_empty_samples{
   my $self = shift;

}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

#sub {
#    my $self = shift;
#}

=head2

 Title   :
 Usage   :
 Function:
 Returns :
 Args    :

=cut

#sub {
#    my $self = shift;
#
#}


=head2

 Title   :
 Usage   :
 Function:
 Returns :
 Args    :

=cut

#sub {
#    my $self = shift;
#
#}

=head2

 Title   :
 Usage   :
 Function:
 Returns :
 Args    :

=cut

#sub {
#    my $self = shift;
#
#}

1;

