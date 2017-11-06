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

#use NGS::Tools::GATK;
use lib '/projects/ccrstaff/jw24/software/TMP/NGS-Tools-GATK/lib/NGS/Tools';
use GATK;

### COMMAND LINE DEFAULT ARGUMENTS ################################################################
our %opts = (
    debug                      => 'N',
    user                       => $ENV{LOGNAME} || $ENV{USER} || getpwuid($<),
    qos                        => 'general-compute',
    partition                  => 'general-compute',
    cluster                    => 'ub-hpc',
	sendemail                  => 'N',
    gatk_version               => '3.8.0',             # set default version
    create_recalibration_table => 'Y',
    create_recalibrated_bam    => 'Y',
    split_bam_by_chr           => 'Y',
    run_haplotype_caller       => 'Y',
    joint_calling_genotypes    => 'N',
    cat_variants               => 'N',
    variant_recalibrator       => 'N',
    apply_recalibration        => 'N',
	merge_hc_bams              => 'N',
    ref       => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/human_g1k_v37_decoy.fasta',
    dbsnp     => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/dbsnp_138.b37.vcf',
    indels    => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/dbsnp_138.b37.vcf',
    hapmap    => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/hapmap_3.3.b37.vcf',
    omni      => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/1000G_omni2.5.b37.vcf',
    g1k_snp   => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/1000G_phase1.snps.high_confidence.b37.vcf',
    g1k_indel => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/1000G_phase1.indels.b37.vcf',
    mills     => '/gpfs/scratch/jw24/GATK/gatk_bundle/2.8/Mills_and_1000G_gold_standard.indels.b37.vcf'
);

GetOptions(
    \%opts,
    "help|?",
    "man",
    "gatk_version|v:s" => \$opts{'gatk_version'},
    "data_type:s",
    "assembly:s",
    "run_id:s",
    "config|c:s",
    "output_dir:s",
    "account:s",
    "user"        => \$opts{'user'},
    "cluster:s"   => \$opts{'cluster'},
    "partition:s" => \$opts{'partition'},
    "qos:s" => \$opts{'qos'},
    "email" => \$opts{'email'},
	"sendemail" => \$opts{'sendemail'},
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
    #"modules:s"                    => \$opts{'modules'},
    #"phaseI"                       => \$opts{'phaseI'},
    #"phaseII"                      => \$opts{'phaseII'},
    #"phaseIII"                     => \$opts{'phaseIII'},
    "create_recalibration_table:s" => \$opts{'create_recalibration_table'},
    "create_recalibrated_bam:s"    => \$opts{'create_recalibrated_bam'},
    "split_bam_by_chr"             => \$opts{'split_bam_by_chr'},
    "run_haplotype_caller:s"       => \$opts{'run_haplotype_caller'},
    "joint_calling_genotypes:s"    => \$opts{'joint_calling_genotypes'},
    "cat_variants:s"               => \$opts{'cat_variants'},
    "variant_recalibrator:s"       => \$opts{'variant_recalibrator'},
    "apply_recalibration:s"        => \$opts{'apply_recalibration'}
) or pod2usage(2);
GetOptions();

####################################################################################################
# Check commandline sanity
####################################################################################################
if (   !$opts{'account'}
    or !$opts{'assembly'}
    or !$opts{'partition'}
    or !$opts{'config'}
    or !$opts{'output_dir'}
    or !$opts{'run_id'}
    or !$opts{'data_type'} )
{
    pod2usage(1);
}

if ( $opts{'help'} ) { pod2usage(1) }
if ( $opts{'man'} ) { pod2usage( -exitstatus => 0, -verbose => 2 ) }

if ( $opts{'assembly'} eq 'GRCh38' ) {

    #    $opts{'ref'} = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/Homo_sapiens_assembly38.fasta';
    $opts{'ref'}       = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna';
    $opts{'dbsnp'}     = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf';
    $opts{'indels'}    = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf';
    $opts{'hapmap'}    = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/hapmap_3.3.hg38.vcf.gz';
    $opts{'omni'}      = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/1000G_omni2.5.hg38.vcf.gz';
    $opts{'g1k_snp'}   = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/1000G_phase1.snps.high_confidence.hg38.vcf.gz';
    $opts{'g1k_indel'} = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz';
    $opts{'mills'}     = '/gpfs/scratch/jw24/GATK/gatk_bundle/GRCh38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz';
}

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
# Preprocessing: base quality recalibration
####################################################################################################

my $base_recalibration = $gatk->run_base_recalibration();

my $recalibrated_bam = $base_recalibration->create_recalibrated_bams();

####################################################################################################
# Variant calling
####################################################################################################

# make variant calls with HaplotypeCaller
# due to long running time and high memory requirement, split recalibrated bams by chr
my $split_bam = $recalibrated_bam->split_bam_by_chr();

my $htc = $split_bam->run_haplotype_caller();

# check if user want to continue for joint calling genotypes
if ( $opts{'joint_calling_genotypes'} eq 'Y' ) {
    my $genotyped_vcf     = $htc->joint_calling_genotypes();
    my $combined_variants = $genotyped_vcf->run_cat_variants();

    ####################################################################################################
    # Variant recalibration
    ####################################################################################################

    my $snp_recalibration              = $combined_variants->run_variant_recalibrator('snp');
    my $apply_variant_recalibrator_snp = $snp_recalibration->run_apply_recalibration('snp');

    my $indel_recalibration               = $apply_variant_recalibrator_snp->run_variant_recalibrator('indel');
    my $apply_variant_recalibration_indel = $indel_recalibration->run_apply_recalibration('indel');

    my $merge_bams = $apply_variant_recalibration_indel->merge_hc_bams();
}
else {
    print "Joint calling genotypes won't be performed\n";
}

# if we reach here, everything is ok
print "All jobs submitted...\n";

exit;

__END__


=head1 NAME

GATK.pl

=head1 DESCRIPTION

B<GATK.pl> Run the Genome Analysis Toolkit pipeline for human samples.


=head1 SYNOPSIS

B<GATK.pl> [options] <--config> <--run_id> <--output_dir> <--account> <--data_type> <--assembly> <--partition>

Arguments:

	--config		A YAML file containing BAM files and sample names for processing
	--run_id		A name/ID for this GATK run, suggest to use the project, dataset or cohort name
	--output_dir		Output directory (for example: /scratch2/user/your_user_name).
	--account		Group name of the user belongs to
	--data_type		Indicate if the input sequence data are from whole genome sequencing or whole exome sequencying. valid values are "wgs", "wes" or "rna-seq"
	--assembly		Assembly version for the ref genome. Use 'b37' or 'GRCh38'. Default to b37
	--partition		Valide values are "general-compute", "industry". Default to "general-compute". 
				If set to "industry", qos and cluster will also be set to "industry"

Options:


	--cluster		Cluster name, valid values are 'ub-hpc', 'industry'. Default to 'ub-hpc'
	--qos			Quality of Server, limits or priority boost. Valid values are "supporters", 
				"general-compute", "industry". Default to "general-compute". 
				Use "industry" if you are submitting jobs to industry cluster
	--email			Used to send Slurm run status. Defaul is "username@buffalo.edu"                       
	--sendemail		Send email after slurm job finishes. Y|N, default is No. Will flood your mailbox if set to 'Y'
	--gatk_version		Version of the GATK program, e.g. "3.3.0", "3.6.0", etc. default "3.8.0"
	--ref			Path for reference genome fasta file (defauld b37)
	--create_recalibration_table		Perform CreateRecalibrationTable (default:Y)
	--merge_raw_vcf		Perfrom MergeRawVcf (default:Y)
	--variant_recalibrator		Perfrom VariantRecalibrator (default:Y)
	--apply_recalibration		Perfrom ApplyRecalibration (default:Y)
	--create_recalibrated_bam		Default is 'Y', yes
	--split_bam_by_chr		Default is 'Y', yes
	--run_haplotype_caller		Default is 'Y', yes
	--joint_calling_genotypes		Default is 'N', no, because users may need to run several batches of samples before making the final genotype calling
	--cat_variants		Defautl is 'N', since default value for --joint_calling_genotypes is set to 'N'. Set to 'Y' if need to perform this step
	--variant_recalibrator		Default is 'N'. Set to 'Y' if need to perform this step
	--apply_recalibration		Default is 'N'. Set to 'Y' if need to perform this step
	--merge_hc_bams		Merge bam files created by HaplotypeCaller for use with IGV. Y|N, Default to No.


Examples:
	
GATK.pl --config config.yaml --output_dir /gpfs/scratch/user/gatk --gatk_version 3.8.0 --data_type wgs --run_id my_project_name --assembly b37 --account your_group name --partition general-compute|industry 
	

=head1 AUTHOR

Jianxin Wang

Center for Computational Research, University at Buffalo

=head1 SEE ALSO

=cut

