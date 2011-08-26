# -*- Perl -*- Test script for sds2expr modules
#

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use lib '.';
use Test::More tests => 8;

BEGIN {
    use_ok('Bio::PCR::Well');
    use_ok('Bio::PCR::Sample');
    use_ok('Bio::PCR::Experiment');
}
require_ok('Bio::PCR::Well');
require_ok('Bio::PCR::Sample');
require_ok('Bio::PCR::Experiment');

## Test sample with real protein data, and ability to add and filter wells
my @wells;
my @ct_vals = (13.566957, 13.569325, 13.875829);

for my $index (0 .. 2) {
    $wells[$index] = Bio::PCR::Well->new(-ct => $ct_vals[$index], -pos => $index+1);
}

my $sample1 = Bio::PCR::Sample->new(-name => "Glut2_1-1000", -wells => \@wells);


@ct_vals = (14.904137, 14.907749, 16.524164);
@wells = ();

for my $index (0 .. 2) {
    $wells[$index] = Bio::PCR::Well->new(-ct => $ct_vals[$index], -pos => $index+1);
}

my $sample2 = Bio::PCR::Sample->new(-name => "Glut2_1-10000", -wells => \@wells);



@ct_vals = (15.61498, 15.55384, 16.557837);
@wells = ();

for my $index (0 .. 2) {
    $wells[$index] = Bio::PCR::Well->new(-ct => $ct_vals[$index], -pos => $index+1);
}

my $sample3 = Bio::PCR::Sample->new(-name => "Glut2_1-20000", -wells => \@wells);


@ct_vals = (15.89332, 16.665163, 16.323917);
@wells = ();

for my $index (0 .. 2) {
    $wells[$index] = Bio::PCR::Well->new(-ct => $ct_vals[$index], -pos => $index+1);
}

my $sample4 = Bio::PCR::Sample->new(-name => "Glut2_1-40000", -wells => \@wells);

my %samples = ($sample1->get_name => $sample1, $sample2->get_name => $sample2, $sample3->get_name => $sample3, $sample4->get_name => $sample4);

my $experiment = Bio::PCR::Experiment->new(-name => "Glut2 Dilutions", -samples => \%samples);

isa_ok($experiment, "Bio::PCR::Experiment", "Is an experiment object");

$experiment->_calc_2ddct_all();

my $exp_with_ref = Bio::PCR::Experiment->new(-name => "Glut2 Dilutions", -samples => \%samples, -ref => $sample4->get_name);

isa_ok($exp_with_ref, "Bio::PCR::Experiment", "Experiment with reference is an experiment object");

$exp_with_ref->_calc_2ddct_all();
