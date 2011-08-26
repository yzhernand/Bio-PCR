# -*- Perl -*- Test script for sds2expr modules
#

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use lib '.';
use Test::More tests => 14;

BEGIN {
    use_ok('Bio::PCR::Well');
    use_ok('Bio::PCR::Sample');
}
require_ok('Bio::PCR::Well');
require_ok('Bio::PCR::Sample');

## Test creation and sample wells
my @ct_vals = (0.1234, 0.2234, 0.3234);
my @wells;

for my $index (0 .. 2) {
    $wells[$index] = Bio::PCR::Well->new(-ct => $ct_vals[$index], -pos => $index+1);
}

my $sample = Bio::PCR::Sample->new(-name => "Test Sample", -wells => \@wells);

isa_ok( $sample, 'Bio::PCR::Sample' );

my @reported_vals = $sample->get_all_ct();

ok( @reported_vals ~~ @ct_vals, "Ct values are correct" );
is( $sample->get_name(), "Test Sample", "Name is correct" );

## Test samples with no name
my $unk_sample = Bio::PCR::Sample->new(-wells => \@wells);
is( $unk_sample->get_name(), "Unknown_sample", "Unknown name is correct" );



## Test sample with real protein data, and ability to add and filter wells
@ct_vals = (13.566957, 13.569325, 13.875829);
@wells = ();

for my $index (0 .. 1) {
    $wells[$index] = Bio::PCR::Well->new(-ct => $ct_vals[$index], -pos => $index+1);
}

my $prot_sample = Bio::PCR::Sample->new(-name => "Glut2");

# Set wells
$prot_sample->set_wells(\@wells);

is($prot_sample->num_wells(), 2, "Sample has 2 wells");

my $new_well = Bio::PCR::Well->new(-ct => $ct_vals[2], -pos => 3);
# Add a well
$prot_sample->add_well($new_well);

is($prot_sample->num_wells(), 3, "Sample now has 3 wells");

@reported_vals = $sample->get_all_ct();

ok( @reported_vals ~~ @ct_vals, "Real data Ct values are correct" );

my $avg = sum(@ct_vals)/3;

is($prot_sample->get_avg_ct(), $avg, "Prefiltered average is correct");

$prot_sample->filter_outliers();

$avg = sum(@ct_vals[0..1])/2;

is($prot_sample->num_wells(), 2, "Filtered sample now has 2 wells");

is($prot_sample->get_avg_ct(), $avg, "Filtered average is correct");

