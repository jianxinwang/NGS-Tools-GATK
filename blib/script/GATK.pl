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

use NGS::Tools::GATK;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
our %opts = (
    debug                      => 'N',
    gatk_version               => '3.8.0',  # set default version
    phaseI                     => 'Y',
    phaseII                    => 'Y',
    phaseIII                   => 'Y',
    create_recalibration_table => 'Y',
    create_recalibrated_bam    => 'Y',
    split_bam_by_chr           => 'Y',
    run_haplotype_caller       => 'Y',
    joint_calling_genotypes    => 'Y',
    cat_variants               => 'Y',
    variant_recalibrator       => 'Y',
    apply_recalibration        => 'Y',
    ref                        => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/human_g1k_v37_decoy.fasta',
    dbsnp                      => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/dbsnp_138.b37.vcf',
    indels                     => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/dbsnp_138.b37.vcf',
    hapmap                     => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/hapmap_3.3.b37.vcf',
    omni                       => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/1000G_omni2.5.b37.vcf',
    g1k_snp                    => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/1000G_phase1.snps.high_confidence.b37.vcf',
    g1k_indel                  => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/1000G_phase1.indels.b37.vcf',
    mills                      => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/Mills_and_1000G_gold_standard.indels.b37.vcf'
#    modules                    => 'samtools,java/1.8.0_45,vcftools,picard/2.2.1,R/3.2.0'
);

GetOptions(
    \%opts,
    "help|?",
    "man",
    "gatk_version|v:s" => \$opts{'gatk_version'},
    "run_id:s",
    "config|c:s",
    "output_dir:s",
    "data_type:s",
    "debug:s"                      => \$opts{'debug'},
    "ref|r:s"                      => \$opts{'ref'},
    "dbsnp:s"                      => \$opts{'dbsnp'},
    "indels:s"                     => \$opts{'indels'},
    "hapmap:s"                     => \$opts{'hapmap'},
    "omni:s"                       => \$opts{'omni'},
    "g1k_snp:s"                    => \$opts{'g1k_snp'},
    "g1k_indel:s"                  => \$opts{'g1k_indel'},
    "mills"                        => \$opts{'mills'},
    "modules:s"                    => \$opts{'modules'},
    "phaseI"                       => \$opts{'phaseI'},
    "phaseII"                      => \$opts{'phaseII'},
    "phaseIII"                     => \$opts{'phaseIII'},
    "create_recalibration_table:s" => \$opts{'create_recalibration_table'},
    "create_recalibrated_bam:s"    => \$opts{'create_recalibrated_bam'},
    "split_bam_by_chr"             => \$opts{'split_bam_by_chr'},
    "run_haplotype_caller:s"       => \$opts{'run_haplotype_caller'},
    "joint_calling_genotypes:s"    => \$opts{'joint_calling_genotypes'},
    "cat_variants:s"               => => \$opts{'cat_variants'},
    "variant_recalibrator:s"       => \$opts{'variant_recalibrator'},
    "apply_recalibration:s"        => \$opts{'apply_recalibration'}
) or pod2usage(2);
GetOptions();

####################################################################################################
# Check commandline sanity
####################################################################################################
if ( !$opts{'config'} or !$opts{'output_dir'} or !$opts{'run_id'} or !$opts{'data_type'}) {
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

####################################################################################################
# Root directories to store output files
####################################################################################################
my $gatk_version  = $opts{'gatk_version'};
my $parent_dir    = File::Spec->rel2abs( $opts{'output_dir'} );
my $phaseI_root   = $parent_dir . '/phaseI_' . $gatk_version;
my $phaseII_root  = $parent_dir . '/phaseII_' . $gatk_version;
my $phaseIII_root = $parent_dir . '/phaseIII_' . $gatk_version;

my $gatk = NGS::Tools::GATK->new(
    bam           => $config,
    config        => $config,
    opts          => \%opts,
    phaseI_root   => $phaseI_root,
    phaseII_root  => $phaseII_root,
    phaseIII_root => $phaseIII_root,
);

####################################################################################################
# Prepare directories for output files
####################################################################################################
# Switch to working directory
chdir($parent_dir);

# Initialize directory structure
$gatk->dir_structure();

foreach my $phasedir ( 'phaseI', 'phaseII', 'phaseIII' ) {
    my $phaseroot = $parent_dir . '/' . $phasedir . '_' . $gatk_version;
    File::Path::make_path( File::Spec->rel2abs($phaseroot) );
    chdir( File::Spec->rel2abs($phaseroot) ) or croak "$!";

    my @phase_dirs = values %{ $gatk->dir->{$phasedir} };
    File::Path::make_path(@phase_dirs);

    # Switch back to parent dir
    chdir( File::Spec->rel2abs($parent_dir) ) or croak "$!";
}

# Fill in the sample_names to the $gatk object
my $sample_names = $gatk->get_sample_name_from_bam();

####################################################################################################
# GATK phase I, indel realigner and base quality recalibration
####################################################################################################

my $base_recalibration = $gatk->run_base_recalibration();

my $recalibrated_bam = $base_recalibration->create_recalibrated_bams();

####################################################################################################
# GATK phase II, variant calling
####################################################################################################

# make variant calls with HaplotypeCaller
# due to long running time and high memory requirement, split recalibrated bams by chr
my $split_bam = $recalibrated_bam->split_bam_by_chr();

my $htc = $split_bam->run_haplotype_caller();

my $genotyped_vcf = $htc->joint_calling_genotypes();

my $combined_variants = $genotyped_vcf->run_cat_variants();

####################################################################################################
# GATK phaseIII, variant recalibration
####################################################################################################

my $snp_recalibration              = $combined_variants->run_variant_recalibrator('snp');
my $apply_variant_recalibrator_snp = $snp_recalibration->run_apply_recalibration('snp');

my $indel_recalibration               = $apply_variant_recalibrator_snp->run_variant_recalibrator('indel');
my $apply_variant_recalibration_indel = $indel_recalibration->run_apply_recalibration('indel');

my $merge_bams = $apply_variant_recalibration_indel->merge_hc_bams();


# if we reach here, everything is ok
print "All jobs submitted...\n";

exit;

__END__


=head1 NAME

GATK.pl

=head1 SYNOPSIS

B<GATK.pl> [options] <--config> <--run_id> <--output_dir>

	Arguments:
	--config				A YAML file containing BAM files and sample names for processing
	--run_id				A name/ID for this GATK run, suggest to use the project, dataset or cohort name
	--output_dir			Output directory (for example: /scratch2/user/your_user_name).
    --data_type             Indicate if the input sequence data are from whole genome sequencing or whole exome sequencying. valid values are "wgs" or "wes"
	
	Options:
	--help					Brief help message
	--man					Full documentation
	--gatk_version			Version of the GATK program, e.g. "3.3.0", "2.4.9", etc. default "3.6.0"
	--ref					Path for reference genome fasta file (default:hg19)
	--indels				Insertion/deletion file in VCF format
	--dbsnp					dbSNP file in VCF format
	--ref_type				Reference genome (default:hg19) Other acceptable reference genome is b37 from the broad institute
	--indel_realigner_target_creator	Perform RealingerTargetCreator (default:Y)
	--indel_realigner			Perform IndelReaglinger (default:Y)
	--create_recalibration_table		Perform CreateRecalibrationTable (default:Y)
	--merge_raw_vcf				Perfrom MergeRawVcf (default:Y)
	--variant_recalibrator			Perfrom VariantRecalibrator (default:Y)
	--apply_recalibration			Perfrom ApplyRecalibration (default:Y)

	Examples:
	
	### Simple case, assuming everying will run smoothly and the default values are fine to you:
	perl gatk.pl --config tumour_normal_paried.yaml --output_dir /scratch2/user/jwang/TEST --gatk_version 3.5.0 --data_type wgs --run_id SIMWGS
	
	
=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--config>

Config YAML file containing samples with corresponding BAM file locations

=item B<--sample>

Name for a job of the GATK pipeline

=item B<--dir>

Output directory name. /scratch/user/ is recommended because this GATK script can generated a few TB of intermediate files if your sample is whole gneome sequencing data.

=item B<--ref>

Path for reference genome fasta file. Default is the path for the b37_decoy reference genome.

=item B<--ref_type>

Reference genome (default:b37_decoy, from Broad Institute, used in 1kg phaseII)

=item B<--pre_processing>

Perform BAM processing including RealignerTargetCreator, IndelRealigner, BaseRecalibrator, PrintReads, and merging and indexing BAM files using Picard.
(default:Y)

=item B<--SNP_calling>

Perform SNP calling including HaplotypeCaller or UnifiedGenotyper, VariantRecalibrator, ApplyRecalibration, and merging VCF files. (default:Y)

=back

=head1 DESCRIPTION

B<GATK.pl> Run the Genome Analysis Toolkit pipeline.

=head1 AUTHOR

Jianxin Wang

Center for Computational Research, University at Buffalo

=head1 SEE ALSO

=cut

