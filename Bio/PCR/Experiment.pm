
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

our $VERSION = '0.5';

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

    my $sample_name = $param{'-name'} // "Unnamed_experiment_";

    my $samples = $param{'-samples'};

    croak "Error when setting samples. '-samples' must be a hash ref\n"
        unless ( ref($samples) eq 'HASH' );

    my $sep   = $param{'-sep'}   // "\t";
    my $ref   = $param{'-ref'}   // undef;
    my $calib = $param{'-calib'} // undef;
    croak "Error: reference sample \"$ref\" not found. ",
        "Please check name is exactly as it appears in input file. Exiting...\n"
        if ( defined($ref) && ( !exists $samples->{$ref} ) );

    $self->{SAMPLE_NAME} = $sample_name;
    $self->{SAMPLES}     = $samples;
    $self->{NUM_SAMPLES} = ( keys %$samples ) * 1;
    $self->{SEP}         = $sep;
    $self->{REF}         = $ref;
    $self->{CALIB}       = $calib;

    bless( $self, $class );
    $self->_discard_empty_samples();

    return $self;
}

=head2 _calc_2ddct

 Title   : _calc_2ddct
 Usage   : $self->_calc_2ddct($ref_sample_name) # INTERNAL METHOD
 Function: Calculates the 2^-ddCt value for the samples given a reference sample name
 Returns : None
 Args    : Reference sample name (string)

=cut

sub _calc_2ddct {
    my ( $self, $reference ) = @_;

    $self->{SAMPLES}->{$reference}->set_delta_ct(0);
    $self->{SAMPLES}->{$reference}->set_dd_ct(0);
    $self->{SAMPLES}->{$reference}->set_2ddct(1);

    for my $target ( keys %{ $self->{SAMPLES} } ) {
        next if ( $target eq $reference );

        my $r_Ct    = $self->{SAMPLES}->{$reference}->get_avg_ct();
        my $t_Ct    = $self->{SAMPLES}->{$target}->get_avg_ct();
        my $t_dct   = $t_Ct - $r_Ct;
        my $t_dd_ct = $t_dct - $self->{SAMPLES}->{$reference}->get_dd_ct();
        my $t_rel_express = 2**( -$t_dd_ct );
        $self->{SAMPLES}->{$target}->set_delta_ct($t_dct);
        $self->{SAMPLES}->{$target}->set_dd_ct($t_dd_ct);
        $self->{SAMPLES}->{$target}->set_2ddct($t_rel_express);
    }

    #say Dumper($self->{SAMPLES});

}

=head2 print_2ddct

 Title   : print_2ddct
 Usage   : $self->print_2ddct();
 Function: Calculate and print the 2^-ddCt value either using a
            reference or using all samples as a reference once
 Returns : None
 Args    : None

=cut

sub print_2ddct {
    my $self = shift;

    my @ref_sample_names
        = ( defined $self->{REF} )
        ? ( $self->{REF} )
        : keys %{ $self->{SAMPLES} };

    say "Results for experiment with sample \"", $self->{SAMPLE_NAME}, "\"";
    say join(
        $self->{SEP},
        qw(Gene_name Sample_name Ct_values Mean_Ct delta_Ct delta_delta_Ct Relative_expression Comment)
    );

    for my $reference (@ref_sample_names) {
        $self->_calc_2ddct($reference);

        my @reference_ct = $self->{SAMPLES}->{$reference}->get_all_ct();

        say join(
            $self->{SEP},
            (   $reference,
                $self->{SAMPLE_NAME},
                $reference_ct[0],
                $self->{SAMPLES}->{$reference}->get_avg_ct(),
                $self->{SAMPLES}->{$reference}->get_delta_ct(),
                $self->{SAMPLES}->{$reference}->get_dd_ct(),
                $self->{SAMPLES}->{$reference}->get_2ddct(),
                "When used as reference"
            )
        );

        for my $ct_pos ( 1 .. @reference_ct - 1 ) {
            say join(
                $self->{SEP},
                ( "", "", $reference_ct[$ct_pos], "", "", "", "", "" )
            );
        }

        for my $target ( keys %{ $self->{SAMPLES} } ) {
            next if ( $target eq $reference );

            my @target_ct = $self->{SAMPLES}->{$target}->get_all_ct();

            say join(
                $self->{SEP},
                (   $target,
                    $self->{SAMPLE_NAME},
                    $target_ct[0],
                    $self->{SAMPLES}->{$target}->get_avg_ct(),
                    $self->{SAMPLES}->{$target}->get_delta_ct(),
                    $self->{SAMPLES}->{$target}->get_dd_ct(),
                    $self->{SAMPLES}->{$target}->get_2ddct(),
                    "When compared to $reference"
                )
            );

            for my $ct_pos ( 1 .. @target_ct - 1 ) {
                say join(
                    $self->{SEP},
                    ( "", "", $target_ct[$ct_pos], "", "", "", "", "" )
                );
            }
        }

        say "";
    }

}

=head2 _discard_empty_samples

 Title   : _discard_empty_samples
 Usage   : $self->_discard_empty_samples # INTERNAL METHOD
 Function: Goes through all samples and discards any which have no methods
 Returns :
 Args    :

=cut

sub _discard_empty_samples {
    my $self = shift;

    for my $sample ( keys $self->{SAMPLES} ) {
        my $num_wells = $self->{SAMPLES}->{$sample}->num_wells;

        $self->{SAMPLES}->{$sample}->filter_outliers()
            if ($num_wells);

        if ( $num_wells == 0 ) {
            delete $self->{SAMPLES}->{$sample};
            $self->{NUM_SAMPLES}--;
        }
    }
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

