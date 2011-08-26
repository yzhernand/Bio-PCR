
=head1 NAME

Bio::PCR::Sample - qPCR Sample object

=head1 SYNOPSIS

    # DO NOT use this module directly. A Sample is part of a
    # Bio::PCR::Experiment object

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

This is a representation of a qPCR sample. It contains
all information on a sample including the wells of the sample,
the average Ct value of all wells, the standard deviation, the
delta Ct, the delta delta Ct, and the 2^(-ddCt) value.

It also includes flags to determine whether or not the sample
is being used as a reference or calibrator, or not.

=head1 FEEDBACK

=head1 AUTHOR - Yozen Hernandez

email - yzhernand at gmail dot com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::PCR::Sample;
use strict;
use 5.010;
use Carp;
use List::Util qw(sum);
use Bio::PCR::Well;

#use Bio::PCR::Experiment;

=head2 new

 Title   : new
 Usage   : my $sample = Bio::PCR::Sample->new(-name => $gene_name, -wells => \@well_list);
 Function: Builds a new Bio::PCR::Sample object
 Returns : Bio::PCR::Sample
 Args    : a hash.  useful keys:
   -name    : The name of the gene/target sequence in this sample
   -wells   : Array reference for a list of Bio::PCR::Well objects.
   -cterr   : Allowed error in Ct value (eg 0.2 allows for an error of +/-0.2 in Ct value)
            When given, this value is used in reference of the sample standard deviation of
            Ct values in the sample.

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = {};

    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param;    # lowercase keys

    my $name = $param{'-name'} // "Unknown_sample";

    #croak "Error: Missing required argument '-wells'."
    #    unless(exists $param{'-wells'});

    my $wells = $param{'-wells'};

    croak "Error: '-wells' not an array ref. ",
         "Please supply an arrayref of Bio::PCR::Wells objects"
        if ( ( defined($wells) ) && ( ref($wells) ne 'ARRAY' ) );

    my $cterr = $param{'-cterr'};

    $self->{NAME}      = $name;
    $self->{WELLS}     = $wells;
    $self->{CT_STDDEV} = $cterr;    # Prefer over sample stddev.

    bless( $self, $class );

    $self->filter_outliers();
    return $self;
}


=head2 get_name

 Title   : get_name
 Usage   : my $gene_name = $sample->get_name()
 Function: Get gene name
 Returns : string, gene name in sample
 Args    : None

=cut

sub get_name {
    my ( $self, $wells ) = @_;

    return $self->{NAME};
}


=head2 set_wells

 Title   : set_wells
 Usage   : $sample->set_wells(\@wells);
 Function: Set the wells in this sample to a specified list
 Returns : Nothing
 Args    : Array ref, containing Bio::PCR::Well objects

=cut

sub set_wells {
    my ( $self, $wells ) = @_;

    ( ref($wells) eq 'ARRAY' )
        ? $self->{WELLS} = $wells
        : croak
        "Error: Using a non-arrayref value when attempting to set the wells.\n";
}


=head2 add_well

 Title   : add_well
 Usage   : $sample->add_well($well);
 Function: Add a well to sample
 Returns : Nothing
 Args    : Bio::PCR::Well objects

=cut

sub add_well {
    my ( $self, $well ) = @_;

    ( ref($well) eq 'Bio::PCR::Well' )
        ? push( @{ $self->{WELLS} }, $well )
        : croak "Error: not a Bio::PCR::Well object when adding a well.\n";
}


=head2 num_wells

 Title   : num_wells
 Usage   : my $num_wells = $sample->num_wells();
 Function: Gets the number of wells in a sample
 Returns : Number of wells
 Args    : None

=cut

sub num_wells {
    my $self = shift;

    return @{ $self->{WELLS} } * 1;
}


=head2 get_all_ct

 Title   : get_all_ct
 Usage   : my @ct_vals = $sample->get_all_ct();
 Function: Returns Ct value for all wells in sample
 Returns : Array of all Ct values
 Args    : None

=cut

sub get_all_ct {
    my $self = shift;

    my @ct_vals;

    for my $well ( @{ $self->{WELLS} } ) {
        push @ct_vals, $well->get_ct();
    }

    return @ct_vals;
}


=head2 get_avg_ct

 Title   : get_avg_ct
 Usage   : my $avg_ct = $sample->get_avg_ct();
 Function: Returns, or calculates and returns, average Ct value for sample
 Returns : Real number, scalar
 Args    : None

=cut

sub get_avg_ct {
    my $self = shift;

    return $self->{AVG_CT}
        if ( defined $self->{AVG_CT} );

    my @ct_vals = $self->get_all_ct();
    my $sum     = sum(@ct_vals);
    $self->{AVG_CT} = $sum / @ct_vals;

    return $self->{AVG_CT};
}


=head2 get_ct_stddev

 Title   : get_ct_stddev
 Usage   : my $stddev_ct = $sample->get_ct_stddev();
 Function: Returns, or calculates and returns, sample
            standard deviation of the Ct values in the sample.
            If the -cterr parameter was given when constructing
            the sample object, that value is returned instead.
 Returns : Real number, scalar
 Args    : None

=cut

sub get_ct_stddev {
    my $self = shift;

    return $self->{CT_STDDEV}
        if ( defined $self->{CT_STDDEV} );

    my $diff_squares = 0;
    for my $val ( $self->get_all_ct() ) {
        $diff_squares += ( $val - $self->{AVG_CT} )**2;
    }

    # Use sample standard deviation (with Bessel's correction)
    # as R, Excel and others do
    my $stddev = sqrt( $diff_squares / ( $self->get_all_ct() - 1 ) );
    $self->{CT_STDDEV} = $stddev;

    return $self->{CT_STDDEV};
}


=head2 filter_outliers

 Title   : filter_outliers
 Usage   : $sample->filter_outliers
 Function: Remove outlying Ct values from Sample, and from consideration
            in average Ct and Ct stddev.
 Returns : None
 Args    : None

=cut

sub filter_outliers() {
    my $self = shift;

    # Calculate average and stddev if they haven't been already
    $self->get_avg_ct
        unless $self->{AVG_CT};
    $self->get_ct_stddev
        unless $self->{CT_STDDEV};

    my @filtered;
    my @wells = @{ $self->{WELLS} };
    for my $well (@wells) {
        my $deviation = abs( $well->get_ct() - $self->{AVG_CT} );
        if ( $deviation <= $self->{CT_STDDEV} ) {
            push @filtered, $well;
        }
        else {
            carp "Warning: value ", $well->get_ct(),
                " lies outside permitted error range. Discarding well at position ",
                $well->get_pos(), " ...\n";

            carp "(|", $well->get_ct(), " - ", $self->{AVG_CT}, "| = ", $deviation , " > $self->{CT_STDDEV})\n";
        }
    }

    croak
        "Error: No values lie within specified error range! Please allow for a larger error! Exiting...\n"
        if ( @filtered == 0 );

    $self->{WELLS} = \@filtered;
    
    # Force recalculation of stats next time they are needed
    $self->{AVG_CT} = undef;
    $self->{CT_STDDEV} = undef;
    $self->get_ct_stddev;
}


=head2 get_delta_ct

 Title   : get_delta_ct
 Usage   : my $delta_ct = $sample->get_delta_ct();
 Function: Returns the delta Ct, if set. Else, returns undef
 Returns : Real number, scalar (if set)
            undef (if unset)
 Args    : None

=cut

sub get_delta_ct() {
    my $self = shift;

    return $self->{D_CT};
}


=head2 set_delta_ct

 Title   : set_delta_ct
 Usage   : $sample->set_delta_ct($delta_ct);
 Function: Sets the delta Ct value for the sample
 Returns : None
 Args    : Real number, the delta Ct value

=cut

sub set_delta_ct() {
    my ( $self, $delta_ct ) = @_;

    # TODO validation
    $self->{D_CT} = $delta_ct;
}


=head2 get_dd_ct

 Title   : get_dd_ct
 Usage   : my $dd_ct = $sample->get_dd_ct();
 Function: Returns the delta delta Ct, if set. Else, returns undef
 Returns : Real number, scalar (if set)
            undef (if unset)
 Args    : None

=cut

sub get_dd_ct() {
    my $self = shift;

    return $self->{DD_CT};
}


=head2 set_dd_ct

 Title   : set_dd_ct
 Usage   : $sample->set_dd_ct($dd_ct);
 Function: Sets the delta delta Ct value for the sample
 Returns : None
 Args    : Real number, the delta delta Ct value

=cut

sub set_dd_ct() {
    my ( $self, $dd_ct ) = @_;

    # TODO validation
    $self->{DD_CT} = $dd_ct;
}


=head2 get_2ddct

 Title   : get_2ddctt
 Usage   : my $rel_express = $sample->get_2ddct();
 Function: Returns the 2^-ddCt (relative expression) value, if set. Else, returns undef.
 Returns : Real number, scalar (if set)
            undef (if unset)
 Args    : None

=cut

sub get_2ddct() {
    my $self = shift;

    return $self->{REL_EXPR};
}


=head2 set_2ddct

 Title   : set_2ddct
 Usage   : $sample->set_2ddct($rel_express);
 Function: Sets the 2^-ddCt value for the sample
 Returns : None
 Args    : Real number, the 2^-ddCt

=cut

sub set_2ddct() {
    my ($self, $rel_express) = @_;

    # TODO validation
    $self->{REL_EXPR} = $rel_express;
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
#sub () {
#    my $self = shift;
#}
#
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
#sub () {
#    my $self = shift;
#}

1;

