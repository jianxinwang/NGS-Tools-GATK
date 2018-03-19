#!/usr/bin/env perl

### GATK.pl ##############################################################################
# Run the Genome Analysis Toolkit pipeline.

### INCLUDES ######################################################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Path::Class;
use File::Path;
use File::Spec;
use File::Basename;
use YAML qw(LoadFile);

#use NGS::Tools::GATK4;
use lib '/user/jw24/GitHub/NGS-Tools-GATK/lib/NGS/Tools';
use GATK4;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
our %opts = ();

GetOptions(
	\%opts,
    "help|?",
    "man",
    "batch_id:s",
    "input_bams:s",
	"config|c:s",
    "output_dir:s",
) or pod2usage(2);

####################################################################################################
# Check commandline sanity
####################################################################################################
 if (   !$opts{'config'}
     or !$opts{'output_dir'}
     or !$opts{'batch_id'}
     or !$opts{'input_bams'} )
 {
     pod2usage(1);
 }
 
if ( $opts{'help'} ) { pod2usage(1) }
if ( $opts{'man'} ) { pod2usage( -exitstatus => 0, -verbose => 2 ) }


# Load the config yaml file
my $config;
if ( $opts{'config'} ) {
    $config = LoadFile( $opts{'config'} );
}
else {
    croak "Config yaml file not found. Please check your commandline\n";
}
# Load the input list yaml file
my $input_bams;
if ( $opts{'input_bams'} ) {
    $input_bams = LoadFile( $opts{'input_bams'} );
}
else {
    croak "Input bam file list not found. Please check your commandline\n";
}

####################################################################################################
# Root directories to store output files
####################################################################################################
my $parent_dir    = File::Spec->rel2abs( $opts{'output_dir'} ); # make output dir optional use hash

my $gatk = NGS::Tools::GATK4->new(
    bam           => $input_bams,
    config        => $config,
	opts          => \%opts,
    phaseI_root   => $parent_dir . '/preprocessing',
    phaseII_root  => $parent_dir . '/variant_discovery',
    phaseIII_root => $parent_dir . '/filtering_annotation',
);

####################################################################################################
# Prepare directories for output files
####################################################################################################
# Switch to working directory
chdir($parent_dir);

# Initialize directory structure
$gatk->dir_structure();

foreach my $phasedir ( $gatk->phaseI_root, $gatk->phaseII_root, $gatk->phaseIII_root ) {
	File::Path::make_path( File::Spec->rel2abs($phasedir) );
    chdir( $phasedir ) or croak "$!";

    my @sub_dirs = values %{ $gatk->subdirs->{$phasedir} };
   	File::Path::make_path(@sub_dirs);
}

# get chromosome names from bam files
$gatk->get_chromosomes();

# get interval list
$gatk->get_intervals();

####################################################################################################
# Preprocessing: base quality recalibration
####################################################################################################
my $base_recalibration = $gatk->run_base_recalibration();
my $recalibrated_bam = $base_recalibration->create_recalibrated_bams();

####################################################################################################
# Variant calling
####################################################################################################
my $htc = $recalibrated_bam->run_haplotype_caller();
my $genomics_db_import = $htc->run_GenomicsDBimport();
my $genotyped_vcf     = $genomics_db_import->joint_calling_genotypes();
my $combined_variants = $genotyped_vcf->gather_vcfs();

####################################################################################################
# Variant recalibration
####################################################################################################
my $snp_recalibration              = $combined_variants->run_variant_recalibrator('SNP');
my $indel_recalibration              = $combined_variants->run_variant_recalibrator('INDEL');
my $apply_variant_recalibrator_snp = $snp_recalibration->run_apply_recalibration('SNP');
my $apply_variant_recalibration_indel = $apply_variant_recalibrator_snp->run_apply_recalibration('INDEL');

# if we reach here, everything is ok
print "All jobs submitted...\n";

exit;

__END__


=head1 NAME

GATK.pl

=head1 DESCRIPTION

B<GATK.pl> Run the Genome Analysis Toolkit pipeline for human samples.


=head1 SYNOPSIS

B<GATK.pl> <--config> <--input_bams> <--batch_id> <--output_dir> 

Arguments:

	--config		A YAML file containing account, cluster, qos information, GATK options and parameters, step to run/skip, etc.
	--input_bams	A YAML file lists input bams
	--batch_id		A identifier for this GATK run
	--output_dir		Output directory (for example: /scratch2/user/your_user_name).


Examples:
	
GATK.pl --config config.yaml --input_bams my_input_bam.yaml --output_dir /gpfs/scratch/user/gatk --batch_id my_first_batch 
	

=head1 AUTHOR

Jianxin Wang

Center for Computational Research, University at Buffalo

=head1 SEE ALSO

=cut

