
=head1 NAME

Bio::PCR::Well - qPCR Well object

=head1 SYNOPSIS

    # DO NOT use this module directly. A well is part of a
    # Bio::PCR::Sample object, which is part of a Bio::PCR::Experiment object

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

This is a description of a qPCR well. It is a fairly simplistic
object, containing only the well position and its Ct value.

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
use warnings;
use Carp;
use 5.010;

#use Bio::PCR::Sample;
#use Bio::PCR::Experiment;

our $VERSION = '0.5';

=head2 new

 Title   : new
 Usage   : my $qpcr_well = Bio::PCR::Well->new(-ct => $ct_val);
 Function: Builds a new Bio::PCR::Well object
 Returns : Bio::PCR::Well
 Args    : a hash.  useful keys:
   -ct      : Ct value for this well
   -pos     : Well position, if any

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = {};

    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param;    # lowercase keys

    # Shouldn't happen: -ct should always be defined.
    my $ct_val = $param{'-ct'} // undef;
    $ct_val = undef
        if ( $ct_val =~ /Undetermined/ );

    #croak "Error: Missing required arugument '-pos'.\n"
    #   unless ($param{'pos');

    my $pos = $param{'-pos'};

    $self->{CT}  = $ct_val;
    $self->{POS} = $pos;

    bless( $self, $class );
    return $self;
}

=head2 get_ct

 Title   : get_ct
 Usage   : my $ct_val = $well->get_ct();
 Function: Retrieve a well's Ct value 
 Returns : Real number, scalar
 Args    : None

=cut

sub get_ct {
    my $self = shift;

    return $self->{CT};
}

=head2 get_pos

 Title   : get_pos
 Usage   : my $pos = $well->get_pos();
 Function: Retrieve a well's position
 Returns : Real number, scalar (if exists)
 Args    : None

=cut

sub get_pos {
    my $self = shift;

    return $self->{POS};
}

1;

