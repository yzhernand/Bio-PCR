# -*- Perl -*- Test script for sds2expr modules
#

use strict;
use warnings;
use Data::Dumper;
use lib '.';
use Test::More tests => 5;

BEGIN {
    use_ok('Bio::PCR::Well');
}
require_ok('Bio::PCR::Well');

my $well = Bio::PCR::Well->new( -ct => 0.5005, -pos => 100 );

isa_ok( $well, 'Bio::PCR::Well' );

is( $well->get_ct(),  0.5005, "Ct value is correct: 0.5005" );
is( $well->get_pos(), 100,    "Position value is correct: 100" );

