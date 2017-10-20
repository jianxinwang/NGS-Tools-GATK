#!perl -T
use 5.006;
use strict;
use warnings;
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'NGS::Tools::GATK' ) || print "Bail out!\n";
}

diag( "Testing NGS::Tools::GATK $NGS::Tools::GATK::VERSION, Perl $], $^X" );
