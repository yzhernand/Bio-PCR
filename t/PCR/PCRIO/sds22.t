# -*- Perl -*- Test script for sds2expr modules
#

use strict;
use Data::Dumper;
use lib '.';
use Test::More tests => 4;

BEGIN {
    use_ok('Bio::PCR::PCRIO::sds22');
}
require_ok('Bio::PCR::PCRIO::sds22');

my $filename = 't/data/test_sds22.txt';

my $sdsfile = Bio::PCR::PCRIO::sds22->new( -file => $filename );

isa_ok( $sdsfile, 'Bio::PCR::PCRIO::sds22' );

my $experiments = $sdsfile->get_all_experiments;

isa_ok( $experiments, 'ARRAY' );

#print Dumper($experiments), "\n";

for my $exp (@$experiments) {
    $exp->print_2ddct();
    print "\n";
}

