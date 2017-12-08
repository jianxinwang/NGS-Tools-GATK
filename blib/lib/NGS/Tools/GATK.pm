package NGS::Tools::GATK;

use strict;
use warnings;
use Moose;
use 5.006;
use Carp;
use File::Basename;
use Params::Validate qw(:all);

use CCR::Utilities::General;
#use lib '/projects/ccrstaff/jw24/software/TMP/CCR-Utilities-General/lib/CCR/Utilities';
#use General2;

=head1 NAME

NGS::Tools::GATK - Genome Analysis Toolkit Perl API

=head1 VERSION

Version 0.05

=cut

our $VERSION = '0.06';

=head1 SYNOPSIS

This package includes a Perl interface to the Genome Analysis Toolkit.

	use NGS::Tools::GATK;

    my $gatk = NGS::Tools::GATK->new(
        bam              => $config,
        ..
    );
    
    ...

=head1 CLASS/INSTANCE VARIABLES

=cut

has 'bam' => ( is => 'ro', );

has 'bam_group' => ( is => 'ro', );

has 'config' => ( is => 'ro', );

=head2 $obj->phaseI_root, $obj->phaseII_root

Root directory for storing phase I and phase II output files

=cut

has 'phaseI_root' => (
    is  => 'ro',
    isa => 'Str'
);

has 'phaseII_root' => (
    is  => 'ro',
    isa => 'Str'
);

has 'phaseIII_root' => (
    is  => 'ro',
    isa => 'Str'
);

=header2

Various SLURM scheduler jobnames

=cut

has 'jobid_base_recalibration' => (
    is     => 'ro',
    writer => 'set_jobid_base_recalibration',
);

has 'jobids_split_gvcf' => (
    is     => 'ro',
    writer => 'set_jobids_split_gvcf',
);

has 'jobid_apply_recalibration' => (
    is => 'rw',

    #isa    => 'Str',
    writer => 'set_jobid_apply_recalibration'
);

has 'jobids_variant_recalibrator' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_jobids_variant_recalibrator'
);

has 'recalibration_table' => (
    is     => 'ro',
    isa    => 'HashRef',
    writer => 'set_recalibration_table',
);

has 'jobids_print_recalibrated_bam' => (
    is     => 'rw',
    writer => 'set_jobids_print_recalibrated_bam'
);

has 'recalibrated_bams' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_recalibrated_bams'
);

has 'sample_names' => (
    is     => 'rw',
    isa    => 'ArrayRef',
    writer => 'set_sample_name',
);

has 'jobids_hc' => (
    is     => 'rw',
    writer => 'set_jobids_hc'
);

has 'output_apply_recalibration' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_output_apply_recalibration',
);

has 'jobid_apply_recalibration' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_jobid_apply_recalibration',
);

has 'recalibrated_vcf' => (
    is => 'rw',

    #isa    => 'Str',
    writer => 'set_recalibrated_vcf'
);

has 'jobid_cat_variants' => (
    is => 'rw',

    #isa    => 'HashRef',
    writer => 'set_jobid_cat_variants'
);

has 'output_cat_variants' => (
    is => 'rw',

    #isa    => 'HashRef',
    writer => 'set_output_cat_variants'
);

=header2

Attributes for output file names

=cut

has 'output_split_bam' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_output_split_bam',
);
has 'jobids_split_bam' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_jobids_split_bam',
);
has 'output_gvcfs' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_output_gvcfs',
);
has 'output_hc_bams' => (
    is     => 'rw',
    isa    => 'HashRef',
    writer => 'set_output_hc_bams',
);

has 'joint_called_genotypes' => (
    is => 'rw',

    #isa    => 'HashRef',
    writer => 'set_joint_called_genotypes'
);

has 'jobid_joint_calling_genotypes' => (
    is => 'rw',

    #isa    => 'HashRef',
    writer => 'set_jobid_joint_calling_genotypes'
);

has 'log_apply_recalibration' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'log/apply_recalibration'
);

has 'output_variant_recalibrator' => (
    is => 'rw',

    #isa    => 'HashRef',
    writer => 'set_output_variant_recalibrator'
);

=head2 $obj->opts

Attributes related to config

=cut

has 'opts' => (
    is  => 'ro',
    isa => 'HashRef',
);

=head2	$obj->covariates
Covariates
=cut

has 'covariates' => (
    is      => 'ro',
    isa     => 'ArrayRef',
    default => sub {
        [qw/ReadGroupCovariate QualityScoreCovariate CycleCovariate ContextCovariate/];
    },
);

=head2 $obj->chromosomes

An array containing chromosomes to be analyzed.

=cut

has 'dir' => (
    is     => 'ro',
    isa    => 'HashRef',
    writer => 'set_dir_structure',
);

=head1 SUBROUTINES/METHODS
=cut

sub dir_structure {
    my $self = shift;

    my %dirs = (
        phaseI => {
            intermediate_bam_files         => 'intermediate_bam_files',
            recalibration_table            => 'recalibration_table',
            scripts                        => 'scripts',
            script_recalibration_table     => 'scripts/recalibration_table',
            script_create_recalibrated_bam => 'scripts/create_recalibrated_bam',
            script_split_bam               => 'scripts/split_bam',
            log                            => 'log',
            log_recalibration_table        => 'log/recalibration_table',
            log_create_recalibrated_bam    => 'log/create_recalibrated_bam',
            log_split_bam                  => 'log/split_bam',
        },
        phaseII => {
            raw_vcf                        => 'raw_vcf',
            scripts                        => 'scripts',
            script_joint_calling_genotypes => 'scripts/joint_calling_genotypes',
            script_haplotype_caller        => 'scripts/haplotype_caller',
            script_cat_variants            => 'scripts/cat_variants',
            script_merge_hc_bams           => 'scripts/merge_hc_bams',
            hc_bam                         => 'hc_bam',
            log                            => 'log',
            log_joint_calling_genotypes    => 'log/joint_calling_genotypes',
            log_haplotype_caller           => 'log/haplotype_caller',
            log_cat_variants               => 'log/cat_variants',
            log_merge_hc_bams              => 'log/merge_hc_bams',
        },
        phaseIII => {
            scripts                     => 'scripts',
            script_variant_recalibrator => 'scripts/variant_recalibrator',
            script_apply_recalibration  => 'scripts/apply_recalibration',
            recalibrated_vcf            => 'recalibrated_vcf',
            vcf_recalibration_model     => 'vcf_recalibration_model',
            log                         => 'log',
            log_variant_recalibrator    => 'log/variant_recalibrator',
            log_apply_recalibration     => 'log/apply_recalibration',
        }
    );

    $self->set_dir_structure( \%dirs );
    $self;
}

=head2 

=cut

sub get_chromosomes {
    my $self = shift;
    my %args = validate(
        @_,
        {
            config => {
                type    => HASHREF,
                default => $self->config()
            }
        }
    );

    my @donors = keys( %{ $args{'config'} } );
    my @bams;    # input bam file with absolute path
    foreach my $d (@donors) {
        if ( exists $args{'config'}->{$d}->{'control'}->{'path'} ) {
            push @bams, $args{'config'}->{$d}->{'control'}->{'path'};
        }
        elsif ( exists $args{'config'}->{$d}->{'case'}->{'path'} ) {
            push @bams, $args{'config'}->{$d}->{'case'}->{'path'};
        }
        else {
            croak "No bam file found with sample \"$d\". Please check your config file";
        }
    }

    # get the header section for just one bam file (assuming all bam files are aligned on the same ref genome)
    my $bam_header = `samtools view -H $bams[0]`;

    my ( %chromosomes, @chromosomes, @length );

    # extract chromosome/contig names from bam header
    while ( $bam_header =~ /\@SQ\s+SN:(.*?)\s+LN:(\d+)/g ) {

        # skip non chromosome contigs
        next if ( $1 =~ /GL|hs37d5|NC_|Un|EBV|random|HLA/i );
        push @chromosomes, $1;
        push @length,      $2;

    }

    for ( my $i = 0 ; $i < ( scalar @chromosomes ) ; $i++ ) {
        $chromosomes{ $chromosomes[$i] } = $length[$i];

        #        last if $i > 3;    # comment out for production
    }

    my %return_value = (
        chromosomes_ordered => \@chromosomes,
        chromosomes         => \%chromosomes,
        donor               => \@donors,
        bam_path            => \@bams
    );

    return ( \%return_value );
}

sub get_gatk_options {
    my $self = shift;
    my %args = validate(
        @_,
        {
            walker => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            step => {
                type     => SCALAR,
                required => 0,
                default  => undef
            }
        }
    );

    my $walker = $args{'walker'};
    my $step   = $args{'step'};

    # initialize with defaul values
    my $time = '72:00:00';
    my $nt   = 1;
    my $nct  = 1;
    my $mem  = 48_000;
    my $rv;
    #my $partition = 'general-compute';

    if ( defined $walker ) {
        if ( $walker eq 'PrintReads' ) {
            $nct  = 1;
            $mem  = 16_000;
            $time = '20:00:00';
        }
        elsif ( $walker eq 'HaplotypeCaller' ) {
            $nct  = 6;
            $mem  = 32_000;
            $time = '72:00:00';
        }
        elsif ( $walker eq 'SelectVariants' ) {
            $nt   = 8;
            $mem  = 48_000;
            $time = '20:00:00';
        }
        elsif ( $walker eq 'GenotypeGVCFs' ) {
            $nt   = 1;
            $mem  = 48_000;
            $time = '72:00:00';
        }
        elsif ( $walker eq 'VariantRecalibrator' ) {
            $time = '20:00:00';
            $nt   = 8;
            $mem  = 48_000;
        }
        elsif ( $walker eq 'ApplyRecalibration' ) {
            $time = '10:00:00';
            $nt   = 8;
            $mem  = 48_000;
        }
        elsif ( $walker eq 'BaseRecalibrator' ) {
            $nct  = 8;
            $mem  = 32_000;
            $time = '10:00:00';
        }
    }
    elsif ( defined $step ) {
        if ( $step eq 'create_recalibrated_bam' ) {
            $time = '30:00:00';
            $nt   = 12;
            $mem  = 12_000;
        }
    }
    $rv->{'nt'}        = $nt;
    $rv->{'nct'}       = $nct;
    $rv->{'mem'}       = $mem;
    $rv->{'time'}      = $time;
    #$rv->{'partition'} = $partition;
    return ($rv);
}

sub get_sample_name_from_bam {
    my $self = shift;
    my @sample_names;

    foreach my $key ( sort keys %{ $self->bam } ) {
        my @bam_files;

        if ( exists $self->bam->{$key}->{'case'}
            and $self->bam->{$key}->{'case'}->{'path'} ne '' )
        {
            push @bam_files, $self->bam->{$key}->{'case'}->{'path'};
        }
        if ( exists $self->bam->{$key}->{'control'}
            and $self->bam->{$key}->{'control'}->{'path'} ne '' )
        {
            push @bam_files, $self->bam->{$key}->{'control'}->{'path'};
        }
        foreach my $bam_file (@bam_files) {
            my $header = `samtools view -H $bam_file`;
            my @header = split( /\n/, $header );
            my %seen;

            foreach my $line (@header) {
                my @read_group_info = split /[\t]/, $line;
                if ( $read_group_info[0] =~ /\@RG/ ) {
                    for ( my $i = 0 ; $i < scalar @read_group_info ; $i++ ) {
                        if ( $read_group_info[$i] =~ /SM:/ ) {    #Obtain SM: tag
                            my ($sample_name) = $read_group_info[$i] =~ /SM:(.*)/;

                            $seen{$sample_name} = ();
                        }
                    }
                }
            }

            #Check if SM names are consistent in the bam header
            if ( scalar keys %seen > 1 ) {
                croak print "sample name in this bam is not consistent. Please fix sample name problem.\n";
            }
            elsif ( scalar keys %seen == 0 ) {
                croak "No sample name found from bam header in \"$bam_file\". Please fix the bam file problem\n";
            }
            push @sample_names, keys %seen;
        }
    }
    $self->set_sample_name( \@sample_names );
    $self;
}

sub make_input_file_string {
    my $self = shift;
    my %args = validate(
        @_,
        {
            input => {
                type     => HASHREF,
                required => 1,
                default  => undef
            },
            sample => {
                type    => SCALAR,
                default => undef,
            }
        }
    );

    my @input_bam_files;

    # get the input bam files
    foreach my $key ( sort keys %{ $args{'input'} } ) {
        push @input_bam_files, $args{'input'}->{$key}->{'control'}->{'path'}
          if ( exists $args{'input'}->{$key}->{'control'}->{'path'} );
        push @input_bam_files, $args{'input'}->{$key}->{'case'}->{'path'}
          if ( exists $args{'input'}->{$key}->{'case'}->{'path'} );
    }

    my $input_bam_files = ' -I ' . join( ' -I ', @input_bam_files );
    return $input_bam_files;
}

sub make_commandline {
    my $self = shift;
    my %args = validate(
        @_,
        {
            memory => {
                type     => SCALAR,
                required => 0,
                default  => 10_000
            },
            walker => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            step => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            input => {

                #type     => SCALAR || hashref,
                required => 0,
                default  => undef
            },
            output => {
                required => 0,
                default  => undef
            },
            tranchesfile => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            recalfile => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            mode => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            covariates => {
                type     => ARRAYREF,
                required => 0,
                default  => undef
            },
            table => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            sample_name => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            nt => {
                type     => SCALAR,
                required => 0,
                default  => 1
            },
            nct => {
                type     => SCALAR,
                required => 0,
                default  => 1
            },
            region => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            bam_output => {
                type     => SCALAR,
                required => 0,
                default  => undef
            },
            data_type => {
                type     => SCALAR,
                required => 0,
                default  => 'wes'
            },

        }
    );

    my $params;
    my $walker       = $args{'walker'};
    my $step         = $args{'step'};
    my $input        = $args{'input'};
    my $output       = $args{'output'};
    my $nt           = $args{'nt'};
    my $nct          = $args{'nct'};
    my $mode         = $args{'mode'};
    my $tranchesfile = $args{'tranchesfile'};
    my $recalfile    = $args{'recalfile'};
    my $covariates   = $args{'covariates'};
    my $total_mem    = $args{'memory'} - 2000;
    my $program      = join( ' ',
        'java', "-Xmx${total_mem}M",
        '-Djava.io.tmpdir=/gpfs/scratch/tmp',
        '-jar $GATK_HOME/GenomeAnalysisTK.jar' );
    my $sample_name = $args{'sample_name'};
    my $table       = $args{'table'};
    my $region      = $args{'region'};
    my $bam_output  = $args{'bam_output'};
    my $data_type   = $args{'data_type'};

    if ( $walker eq 'BaseRecalibrator' ) {
        my $input_bam_files = $input;

        my $known_sites = ' -knownSites '
          . join( ' -knownSites ',
            $self->opts->{'mills'},  $self->opts->{'g1k_indel'}, $self->opts->{'dbsnp'},
            $self->opts->{'hapmap'}, $self->opts->{'g1k_snp'} );

        my $input_covariates = '';
        foreach my $covariate ( @{$covariates} ) {
            $input_covariates = join( ' ', $input_covariates, "--covariate", $covariate );
        }

        $params = join( ' ',
            $input_bam_files, '-R', $self->opts->{'ref'},
            '-l INFO', '-nct', $nct, '-nt', $nt, $known_sites, $input_covariates, '--bqsrBAQGapOpenPenalty', 30, '-o',
            $output );

    }
    elsif ( $walker eq 'PrintReads' ) {
        if ( $step eq 'create_recalibrated_bam' ) {
            $params = join( ' ',
                $input, '-R', $self->opts->{'ref'},
                '-BQSR', $table, '-nct', $nct, '--disable_indel_quals', '-o', $output );

        }
        elsif ( $step eq 'split_bam_by_chr' ) {
            $params = join( ' ', $input, '-R', $self->opts->{'ref'}, '-L', $region, '-nct', $nct, '-o', $output );
        }
    }
    elsif ( $walker eq 'HaplotypeCaller' ) {
		
        $params = join( ' ',
            '-R', $self->opts->{'ref'},
            '-I',          $input, '-o',   $output, '-l', 'INFO', "--allow_potentially_misencoded_quality_scores",
            '-minPruning', 2,      '-nct', 1,       '--dbsnp',
            $self->opts->{'dbsnp'}, '--bamOutput', $bam_output,
            '--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000' );
    }
    elsif ( $walker eq 'GenotypeGVCFs' ) {

        $params = join( ' ',
            '-R', $self->opts->{'ref'},
            '-V', $input, '-o', $output, '-newQual', '-l', 'INFO', '--dbsnp', $self->opts->{'dbsnp'},
            '-L', $region );

    }
    elsif ( $walker eq 'VariantRecalibrator' ) {
        my $resource;

        if ( $mode eq 'snp' ) {
            $resource = join( ' ',
                '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', $self->opts->{'hapmap'},
                '-resource:omni,known=false,training=true,truth=true,prior=12.0',   $self->opts->{'omni'},
                '-resource:1000G,known=false,training=true,truth=false,prior=10.0', $self->opts->{'g1k_snp'},
                '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0',  $self->opts->{'dbsnp'},
            );
        }
        elsif ( $mode eq 'indel' ) {
            $resource = join( ' ',
                '-resource:mills,known=false,training=true,truth=true,prior=12.0', $self->opts->{'mills'},
                '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0', $self->opts->{'dbsnp'} );
        }

        my $params_common = join( ' ',
            '-R',
            $self->opts->{'ref'},
            '-nt',
            $nt,
            '-an QD',
            '-an FS',
            '-an SOR',
            '-an MQ',
            '-an MQRankSum',
            '-an ReadPosRankSum',
            '-an InbreedingCoeff',
            '-tranche 100.0',
            '-tranche 99.95',
            '-tranche 99.9',
            '-tranche 99.5',
            '-tranche 99.0',
            '-tranche 97.0',
            '-tranche 96.0',
            '-tranche 95.0',
            '-tranche 94.0',
            '-tranche 93.0',
            '-tranche 92.0',
            '-tranche 91.0',
            '-tranche 90.0',
            '-input',
            $input,
            $resource,
            '-recalFile',
            $output->{$mode}->{'recal'},
            '-tranchesFile',
            $output->{$mode}->{'tranches'},
            '-rscriptFile',
            $output->{$mode}->{'rscript'} );

        if ( $data_type eq 'wgs' ) {
            $params_common =~ s/\-an QD/\-an QD \-an DP/;
        }
        elsif ( $data_type eq 'wes' ) {
            $params_common =~ s/\-an InbreedingCoeff//;
			
        }
		elsif ( $data_type eq 'rna-seq' ) {
			$params_common =~ s/\-an InbreedingCoeff/-dontUseSoftClippedBases -stand_call_conf 20.0/;
		}

        if ( $mode eq 'snp' ) {
            $params = join( ' ', $params_common, '-mode SNP', '--maxGaussians 6' );
        }
        elsif ( $mode eq 'indel' ) {
            $params = join( ' ', $params_common, '-mode INDEL', '--maxGaussians 6' );
        }
    }
    elsif ( $walker eq 'ApplyRecalibration' ) {
        my $params_common = join( ' ',
            '-R',          $self->opts->{'ref'}, '-input',   $input, '-tranchesFile',
            $tranchesfile, '-recalFile',         $recalfile, '-o',   $output );

        if ( $mode eq 'snp' ) {
            $params = join( ' ', $params_common, '-mode SNP', '--ts_filter_level 99.6' );
        }
        elsif ( $mode eq 'indel' ) {
            $params = join( ' ', $params_common, '-mode INDEL', '--ts_filter_level 95.0' );
        }

    }
    elsif ( $walker eq 'SelectVariants' ) {

        #my $sample_arg = join( ' ', '-sn', $sample_name );
        $params = join( ' ', '-R', $self->opts->{'ref'}, '-V', $input, '-o', $output, '-L', $region );
    }

    my $cmd = join( ' ', $program, "-T $walker", $params );
    my $check_file = "while [ ! -f \$GATK_HOME/GenomeAnalysisTK.jar ]
do sleep 5
done
    ";

    #$cmd = "\n" . $check_file . "\n" . $cmd;
    print "$cmd\n";
    return ($cmd);

}

sub run_base_recalibration {
    my $self = shift;
    my $jobids_base_recalib;
    my $recal_tables;

    print "\nCreating \'Base Quality Recalibration Table\' command...\n";

    my $gatk_options = $self->get_gatk_options( walker => 'BaseRecalibrator' );

    foreach my $donor ( sort keys %{ $self->config } ) {
		foreach my $sample_type (qw(case control)) {
	        my $sample_id   = $self->config->{$donor}->{$sample_type}->{'ID'};
			next if (! $sample_id); # the specific $smaple_type does exist, so skip
	        my $input_bam   = $self->config->{$donor}->{$sample_type}->{'path'};
	        my $output_file = join( '_', $self->opts->{'run_id'}, $sample_id, 'recalibration_table.grp' );
	        my $output_dir  = join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'recalibration_table'} );
	        $output_file = join( '/', $output_dir, $output_file );
	
	        my $commandline = $self->make_commandline(
	            walker     => 'BaseRecalibrator',
	            input      => "-I " . $input_bam,
	            covariates => $self->covariates,
	            nct        => $gatk_options->{'nct'},
	            nt         => 1,
	            memory     => $gatk_options->{'mem'},
	            output     => $output_file,
	        );
	
	        my $jobname_base_recalibration = join( '_', 'BQR', $sample_id, $self->opts->{'run_id'} );
	        my $scriptname_base_recalibration = 'slurm_' . $jobname_base_recalibration . ".sh";
	        my $slurm_jobid_crt;

	        if ( $self->opts->{'create_recalibration_table'} eq 'Y' ) {
	            $slurm_jobid_crt = CCR::Utilities::General::create_submit_slurm_job(
	                debug      => $self->opts->{'debug'},
	                account    => $self->opts->{'account'},
	                email      => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
	                sendemail  => $self->opts->{'sendemail'},
	                partition  => $self->opts->{'partition'},
	                qos        => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
	                cluster    => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
	                script     => $scriptname_base_recalibration,
	                command    => $commandline,
	                script_dir => join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'script_recalibration_table'} ),
	                modules    => "java/1.8.0_45,gatk/" . $self->opts->{'gatk_version'},
	                job_name   => $jobname_base_recalibration,
	                time       => $gatk_options->{'time'},
	                nodes      => 1,
	                memory     => $gatk_options->{'mem'},
	                ntasks_per_node => $gatk_options->{'nct'},
	                dependency      => 'none',
	                logfile          => join( '/',
	                    $self->phaseI_root, $self->dir->{'phaseI'}->{'log_recalibration_table'},
	                    $jobname_base_recalibration )
	            );
	
	            $jobids_base_recalib->{$sample_id} = $slurm_jobid_crt;
	        }
	        else {
	            print "BaseQualRecalibration will not be performed!!\n";
	        }
	
	        $recal_tables->{$sample_id} = $output_file;
	    	$self->set_jobid_base_recalibration($jobids_base_recalib) if ($jobids_base_recalib);
	    	$self->set_recalibration_table($recal_tables);
		}
	}
	$self;
}

sub create_recalibrated_bams {
    my $self = shift;
    my $jobids_create_recalibrated_bam;
    my $jobname_create_recalibrated_bam;
    my $recalibrated_bams;
    my $gatk_options = $self->get_gatk_options( step => 'create_recalibrated_bam' );

    foreach my $donor ( sort keys %{ $self->config } ) {
		foreach my $sample_type (qw(case control) ) {
	        my $sample_id = $self->config->{$donor}->{$sample_type}->{'ID'};
			next if (! $sample_id);
	        my $input_bam = $self->config->{$donor}->{$sample_type}->{'path'};
	
	        my $dir = join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'intermediate_bam_files'} );
	        my $output_file = join( '_', $self->opts->{'run_id'}, $sample_id, 'recalibrated.bam' );
	        $output_file = join( '/', $dir, $output_file );
	
	        # make commandline for slurm script
	        my $commandline = $self->make_commandline(
	            walker      => 'PrintReads',
	            step        => 'create_recalibrated_bam',
	            input       => "-I " . $input_bam,
	            output      => $output_file,
	            table       => $self->recalibration_table->{$sample_id},
	            sample_name => $sample_id,
	            nct         => $gatk_options->{'nct'},
	            memory      => $gatk_options->{'mem'}
	        );
	
	        my $scriptname_create_recalibrated_bam =
	          join( '_', 'slurm', 'BaseQualityRecalibratedBam', $sample_id, $self->opts->{'run_id'} . "\.sh" );
	        my $jobname_create_recalibrated_bam = join( '_', 'BQRB', $sample_id, $self->opts->{'run_id'} );
	
	        my $slurm_jobid_create_recalibrated_bam;
	        my $dependency = 'none';
	        if ( $self->opts->{'create_recalibration_table'} eq 'Y' ) {
	            $dependency = $self->jobid_base_recalibration->{$sample_id};
	        }
	
	        if ( $self->opts->{'create_recalibrated_bam'} eq 'Y' ) {
	            $slurm_jobid_create_recalibrated_bam = CCR::Utilities::General::create_submit_slurm_job(
	                debug     => $self->opts->{'debug'},
	                account   => $self->opts->{'account'},
	                email     => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
	                sendemail => $self->opts->{'sendemail'},
	                partition => $self->opts->{'partition'},
	                qos       => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
	                cluster   => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
	                script    => $scriptname_create_recalibrated_bam,
	                command   => $commandline,
	                script_dir =>
	                  join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'script_create_recalibrated_bam'} ),
	                modules         => "java/1.8.0_45,gatk/" . $self->opts->{'gatk_version'},
	                job_name        => $jobname_create_recalibrated_bam,
	                time            => $gatk_options->{'time'},
	                nodes           => 1,
	                memory          => $gatk_options->{'mem'},                                  #12000,
	                ntasks_per_node => $gatk_options->{'nct'},                                  #4,
	                dependency      => $dependency,
	                logfile          => join( '/',
	                    $self->phaseI_root, $self->dir->{'phaseI'}->{'log_create_recalibrated_bam'},
	                    $jobname_create_recalibrated_bam )
	            );
	        }
	        else {
	            print "PrintReads will not be performed!!\n";
	        }
	
	        $recalibrated_bams->{$sample_id} = $output_file;
	
	        $jobids_create_recalibrated_bam->{$sample_id} = $slurm_jobid_create_recalibrated_bam
	          if ( defined $slurm_jobid_create_recalibrated_bam );
	
	    	$self->set_recalibrated_bams($recalibrated_bams);
	    	$self->set_jobids_print_recalibrated_bam($jobids_create_recalibrated_bam);
		}
	}
    $self;
}

sub split_bam_by_chr {
    my $self         = shift;
    my $gatk_options = $self->get_gatk_options( walker => 'PrintReads' );
    my $output_dir   = join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'intermediate_bam_files'} );
    my $jobid_split_bam_;
    my $output_file_split_bam_;
    my $chromosome = $self->get_chromosomes->{'chromosomes'};

    foreach my $donor ( sort keys %{ $self->config } ) {
		foreach my $sample_type (qw(case control) ) {
	        my $sample_id = $self->config->{$donor}->{$sample_type}->{'ID'};
			next if (! $sample_id);
	        foreach my $chr ( sort keys %$chromosome ) {
	            my $output_file =
	              join( '/', $output_dir, $self->opts->{'run_id'} . '_' . $sample_id . '_' . $chr . '.bam' );
	            my $input       = $self->recalibrated_bams->{$sample_id};
	            my $commandline = $self->make_commandline(
	                walker      => 'PrintReads',
	                step        => 'split_bam_by_chr',
	                input       => '-I ' . $input,
	                output      => $output_file,
	                nt          => $gatk_options->{'nt'},
	                memory      => $gatk_options->{'mem'},
	                sample_name => $sample_id,
	                region      => $chr,
	            );
	            my $jobname_split_bam    = 'SplitBam_' . $sample_id . '_' . $chr . '_' . $self->opts->{'run_id'};
	            my $scriptname_split_bam = 'slurm_' . $jobname_split_bam . '.sh';
	            my $dependency           = 'none';
	            if ( $self->opts->{'create_recalibrated_bam'} eq 'Y' ) {
	                $dependency = $self->jobids_print_recalibrated_bam->{$sample_id};
	            }
	
	            my $slurm_jobid_split_bam;
	
	            if ( $self->opts->{'split_bam_by_chr'} eq 'Y' ) {
	                $slurm_jobid_split_bam = CCR::Utilities::General::create_submit_slurm_job(
	                    debug           => $self->opts->{'debug'},
	                    account         => $self->opts->{'account'},
	                    email           => $self->opts->{'email'}|| $self->opts->{'user'}. '@buffalo.edu',
	                    sendemail       => $self->opts->{'sendemail'},
	                    partition       => $self->opts->{'partition'},
	                    qos             => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
	                    cluster         => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
	                    script          => $scriptname_split_bam,
	                    command         => $commandline,
	                    script_dir      => join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'script_split_bam'} ),
	                    modules         => "java/1.8.0_45,gatk/" . $self->opts->{'gatk_version'},
	                    job_name        => $jobname_split_bam,
	                    time            => $gatk_options->{'time'},
	                    nodes           => 1,
	                    memory          => $gatk_options->{'mem'},
	                    ntasks_per_node => $gatk_options->{'nct'},
	                    dependency      => $dependency,
	                    logfile         =>
	                      join( '/', $self->phaseI_root, $self->dir->{'phaseI'}->{'log_split_bam'}, $jobname_split_bam )
	                );
	            }
	            $jobid_split_bam_->{$sample_id}->{$chr}       = $slurm_jobid_split_bam;
	            $output_file_split_bam_->{$sample_id}->{$chr} = $output_file;
	        }
	   
	    	$self->set_output_split_bam($output_file_split_bam_);
	    	$self->set_jobids_split_bam($jobid_split_bam_);
		}
	}
    $self;
}

sub run_haplotype_caller {
    my $self = shift;
    my $jobids_haplotype_caller;
    my $output_gvcfs;
    my $output_hc_bams;
    my $output_dir = join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'raw_vcf'} );
    my $gatk_options = $self->get_gatk_options( walker => 'HaplotypeCaller' );
    my $chromosome = $self->get_chromosomes->{'chromosomes'};

    foreach my $chr ( sort keys %$chromosome ) {
        foreach my $donor ( sort keys %{ $self->config } ) {
			foreach my $sample_type (qw(case control) ) {
	            my $sample_id      = $self->config->{$donor}->{$sample_type}->{'ID'};
				next if (! $sample_id);
	            my $bam_file       = $self->output_split_bam->{$sample_id}->{$chr};
	            my $output_file    = join( '/', $output_dir, basename($bam_file) . '.g.vcf' );
	            my $bam_outputfile = $output_file;
	            $bam_outputfile =~ s/\.bam\.g\.vcf/_HC\.bam/;
	            $bam_outputfile =~ s/raw_vcf/hc_bam/;
	
	            my $commandline = $self->make_commandline(
	                walker     => 'HaplotypeCaller',
	                input      => $bam_file,
	                output     => $output_file,
	                bam_output => $bam_outputfile,
	                memory     => $gatk_options->{'mem'},
	                nct        => $gatk_options->{'nct'},
	            );
	
	            my $jobname_haplotype_caller = join( '_', 'HC', $sample_id, $chr, $self->opts->{'run_id'} );
	            my $scriptname_haplotype_caller = 'slurm_' . $jobname_haplotype_caller . "\.sh";
	
	            my $dependency = 'none';
	            if ( $self->opts->{'split_bam_by_chr'} eq 'Y' ) {
	                $dependency = $self->jobids_split_bam->{$sample_id}->{$chr};
	            }
	
	            my $slurm_jobid_hc;
	            if ( $self->opts->{'run_haplotype_caller'} eq 'Y' ) {
	                $slurm_jobid_hc = CCR::Utilities::General::create_submit_slurm_job(
	                    debug			=> $self->opts->{'debug'},
	             		account         => $self->opts->{'account'},
	                    email           => $self->opts->{'email'}|| $self->opts->{'user'}. '@buffalo.edu',
	                    sendemail       => $self->opts->{'sendemail'},
	                    partition       => $self->opts->{'partition'},
	                    qos             => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
	                    cluster         => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
	                    script  => $scriptname_haplotype_caller,
	                    command => $commandline,
	                    script_dir =>
	                      join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'script_haplotype_caller'} ),
	                    modules         => "java/1.8.0_45,gatk/" . $self->opts->{'gatk_version'},
	                    job_name        => $jobname_haplotype_caller,
	                    time            => $gatk_options->{'time'},
	                    nodes           => 1,
	                    memory          => $gatk_options->{'mem'},
	                    ntasks_per_node => $gatk_options->{'nct'},
	                    dependency      => $dependency,
	                    logfile          => join( '/',
	                        $self->phaseII_root, $self->dir->{'phaseII'}->{'log_haplotype_caller'},
	                        $jobname_haplotype_caller )
	                );
	            }
	            $jobids_haplotype_caller->{$sample_id}->{$chr} = $slurm_jobid_hc;
	            $output_gvcfs->{$sample_id}->{$chr}            = $output_file;
	            $output_hc_bams->{$sample_id}->{$chr}          = $bam_outputfile;
	        }    # end of foreach $bam file
	    
	    	$self->set_output_gvcfs($output_gvcfs);
	    	$self->set_output_hc_bams($output_hc_bams);
	    	$self->set_jobids_hc($jobids_haplotype_caller);
		}
	}
    $self;
}

sub joint_calling_genotypes {
    my $self = shift;
    my $joint_called_genotypes;
    my $jobid_joint_calling_genotypes;

    my $chromosome   = $self->get_chromosomes->{'chromosomes'};
    my $gatk_options = $self->get_gatk_options( walker => 'GenotypeGVCFs' );
    my $output_dir   = join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'raw_vcf'} );
    foreach my $chr ( sort keys %$chromosome ) {
        my @gvcfs;
        my @dependency;
        foreach my $donor ( sort keys %{ $self->config } ) {
			foreach my $sample_type (qw(case control) ) {
            	my $sample_id = $self->config->{$donor}->{$sample_type}->{'ID'};
				next if (! $sample_id);
            	push @gvcfs,      $self->output_gvcfs->{$sample_id}->{$chr};
            	push @dependency, $self->jobids_hc->{$sample_id}->{$chr};
        	}
		}
        my $output_file = join( '/', $output_dir, $self->opts->{'run_id'} . '_' . $chr . '_raw.vcf' );
        my $gvcfs = join( " -V ", @gvcfs );

        my $commandline = $self->make_commandline(
            walker => 'GenotypeGVCFs',
            input  => $gvcfs,
            output => $output_file,
            nt     => $gatk_options->{'nt'},
            memory => $gatk_options->{'mem'},
            region => $chr,
        );

        my $jobname_joint_calling_genotypes    = 'GenotypeGVCF_' . $chr . '_' . $self->opts->{'run_id'};
        my $scriptname_joint_calling_genotypes = 'slurm_' . $jobname_joint_calling_genotypes . '.sh';

        my $dependency                          = 'none';
        my $slurm_jobid_joint_calling_genotypes = 'none';

        if ( $self->opts->{'run_haplotype_caller'} eq 'Y' ) {
            $dependency = join( ':', @dependency );
        }
        if ( $self->opts->{'joint_calling_genotypes'} eq 'Y' ) {
            $slurm_jobid_joint_calling_genotypes = CCR::Utilities::General::create_submit_slurm_job(
                debug     => $self->opts->{'debug'},
                account   => $self->opts->{'account'},
                email     => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
                sendemail => $self->opts->{'sendemail'},
                partition => $self->opts->{'partition'},
                qos       => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
                cluster   => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
                script    => $scriptname_joint_calling_genotypes,
                command   => $commandline,
                script_dir =>
                  join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'script_joint_calling_genotypes'} ),
                modules         => "java/1.8.0_45,gatk/" . $self->opts->{'gatk_version'},
                job_name        => $jobname_joint_calling_genotypes,
                time            => $gatk_options->{'time'},
                nodes           => 1,
                memory          => $gatk_options->{'mem'},
                ntasks_per_node => $gatk_options->{'nct'},
                dependency      => $dependency,
                logfile          => join( '/',
                    $self->phaseII_root, $self->dir->{'phaseII'}->{'log_joint_calling_genotypes'},
                    $jobname_joint_calling_genotypes )
            );
        }
        $joint_called_genotypes->{$chr}        = $output_file;
        $jobid_joint_calling_genotypes->{$chr} = $slurm_jobid_joint_calling_genotypes;
    }
    $self->set_jobid_joint_calling_genotypes($jobid_joint_calling_genotypes);
    $self->set_joint_called_genotypes($joint_called_genotypes);
    $self;
}

sub run_cat_variants {
    my $self         = shift;
    my $gatk_options = $self->get_gatk_options( walker => 'CatVariants', );
    my $chromosome   = $self->get_chromosomes->{'chromosomes_ordered'};
    my @input_vcfs;
    my @depends;
    foreach my $chr ( @{$chromosome} ) {
        push @input_vcfs, $self->joint_called_genotypes->{$chr};
        push @depends,    $self->jobid_joint_calling_genotypes->{$chr};
    }
    my $input_vcf   = join( " -V ", @input_vcfs );
    my $output_dir  = join( '/',    $self->phaseII_root, $self->dir->{'phaseII'}->{'raw_vcf'} );
    my $output_file = join( '/',    $output_dir, $self->opts->{'run_id'} . '_merged_raw.vcf' );
    my $commandline = join( " ",
        "java -Xmx5000m -cp \$GATK_HOME/GenomeAnalysisTK.jar",
        "org.broadinstitute.gatk.tools.CatVariants",
        "-R", $self->opts->{'ref'},
        "-V", $input_vcf, "-out", $output_file, "-assumeSorted" );
    my $jobname_cat_variants     = 'CatVariant' . $self->opts->{'run_id'};
    my $scriptname_cat_variants  = 'slurm_' . $jobname_cat_variants . '.sh';
    my $dependency               = join( ':', @depends );
    my $slurm_jobid_cat_variants = 'none';

    if ( $self->opts->{'cat_variants'} eq 'Y' ) {

        $slurm_jobid_cat_variants = CCR::Utilities::General::create_submit_slurm_job(
            debug           => $self->opts->{'debug'},
            account         => $self->opts->{'account'},
            email           => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
            sendemail       => $self->opts->{'sendemail'},
            partition       => $self->opts->{'partition'},
            qos             => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
            cluster         => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
            script          => $scriptname_cat_variants,
            command         => $commandline,
            script_dir      => join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'script_cat_variants'} ),
            modules         => "java/1.8.0_45,gatk/" . $self->opts->{'gatk_version'},
            job_name        => $jobname_cat_variants,
            time            => $gatk_options->{'time'},
            nodes           => 1,
            memory          => $gatk_options->{'mem'},
            ntasks_per_node => $gatk_options->{'nt'},
            dependency      => $dependency,
            logfile         =>
              join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'log_cat_variants'}, $jobname_cat_variants )
        );
    }
    $self->set_jobid_cat_variants($slurm_jobid_cat_variants);
    $self->set_output_cat_variants($output_file);
    $self;
}

sub run_variant_recalibrator {
    my ( $self, $mode ) = @_;
    my $output_file;
    my $jobids_variant_recalibrator;
    my $gatk_options = $self->get_gatk_options( walker => 'VariantRecalibrator', );

    print "\nBuilding \'recalibration\' model for $mode ... \n";
    my $input_file = $self->output_cat_variants;
    if ( $mode eq 'indel' ) {
        $input_file = $self->recalibrated_vcf->{'snp'};
    }
    my $output_dir = join( '/', $self->phaseIII_root, $self->dir->{'phaseIII'}->{'vcf_recalibration_model'} );
    my $recal_file = join( '/', $output_dir,          $self->opts->{'run_id'} . '_output_' . $mode . '.recal' );

    my $tranches_file = join( '/', $output_dir, $self->opts->{'run_id'} . '_output_' . $mode . '.tranches' );
    my $rscript_file  = join( '/', $output_dir, $self->opts->{'run_id'} . '_output_' . $mode . '.plots.R' );

    $output_file->{$mode}->{'tranches'} = $tranches_file;
    $output_file->{$mode}->{'recal'}    = $recal_file;
    $output_file->{$mode}->{'rscript'}  = $rscript_file;

    my $commandline = $self->make_commandline(
        walker    => 'VariantRecalibrator',
        mode      => $mode,
        data_type => $self->opts->{'data_type'},
        input     => $input_file,
        output    => $output_file,
        nt        => $gatk_options->{'nt'},
        memory    => $gatk_options->{'mem'},
    );

    my $jobname_variant_recalibrator = join( '_', 'VRC', $mode, $self->opts->{'run_id'} );
    my $scriptname_variant_recalibrator = 'slurm_' . $jobname_variant_recalibrator . '.sh';

    my $slurm_jobid_variant_recalibrator = 'none';
    my $dependency                       = 'none';

    if ( $mode eq 'snp' ) {
        $dependency = $self->jobid_cat_variants;
    }
    elsif ( $mode eq 'indel' ) {
        $dependency = $self->jobid_apply_recalibration->{'snp'};
    }

    if ( $self->opts->{'variant_recalibrator'} eq 'Y' ) {
        $slurm_jobid_variant_recalibrator = CCR::Utilities::General::create_submit_slurm_job(
            debug      => $self->opts->{'debug'},
            account    => $self->opts->{'account'},
            email      => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
            sendemail  => $self->opts->{'sendemail'},
            partition  => $self->opts->{'partition'},
            qos        => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
            cluster    => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
            script     => $scriptname_variant_recalibrator,
            command    => $commandline,
            script_dir => join( '/', $self->phaseIII_root, $self->dir->{'phaseIII'}->{'script_variant_recalibrator'} ),
            modules    => "java/1.8.0_45,R,gatk/" . $self->opts->{'gatk_version'},
            job_name   => $jobname_variant_recalibrator,
            time       => $gatk_options->{'time'},
            nodes      => 1,
            memory     => $gatk_options->{'mem'},
            ntasks_per_node => $gatk_options->{'nt'},
            dependency      => $dependency,
            logfile          => join( '/',
                $self->phaseIII_root, $self->dir->{'phaseIII'}->{'log_variant_recalibrator'},
                $jobname_variant_recalibrator )
        );

    }
    else {
        print "VariantRecalibrator will not be performed!!\n";
    }

    $jobids_variant_recalibrator->{$mode} = $slurm_jobid_variant_recalibrator;
    $self->set_jobids_variant_recalibrator($jobids_variant_recalibrator);
    $self->set_output_variant_recalibrator($output_file);
    $self;
}

sub run_apply_recalibration {
    my ( $self, $mode ) = @_;
    my $gatk_options = $self->get_gatk_options( walker => 'ApplyRecalibration', );
    my $output_apply_recalibration;
    my $jobid_apply_recalibration;

    my $jobnames_apply_recalibration;
    my $recalibrated_vcf;
    my $slurm_jobid_apply_recalibration = 'none';
    my $input_vcf                       = $self->output_cat_variants;    #$self->joint_called_genotypes->{$chr};
    if ( $mode eq 'indel' ) {
        $input_vcf = $self->recalibrated_vcf->{'snp'};
    }

    my $output_dir = join( '/', $self->phaseIII_root, $self->dir->{'phaseIII'}->{'recalibrated_vcf'} );
    my $output_file = $input_vcf;
    $output_file =~ s/\.vcf/_${mode}_recalibrated\.vcf/;
    $output_file = join( '/', $output_dir, basename($output_file) );

    my $commandline = $self->make_commandline(
        walker       => 'ApplyRecalibration',
        mode         => $mode,
        input        => $input_vcf,
        output       => $output_file,
        tranchesfile => $self->output_variant_recalibrator->{$mode}->{'tranches'},
        recalfile    => $self->output_variant_recalibrator->{$mode}->{'recal'},
        memory       => $gatk_options->{'mem'},
    );
    $jobnames_apply_recalibration = join( '_', 'ARCal', $mode, $self->opts->{'run_id'} );
    my $scriptname_apply_recalibration = 'slurm_' . $jobnames_apply_recalibration . '.sh';

    my $dependency = 'none';

    if ( $self->jobids_variant_recalibrator ) {
        $dependency = $self->jobids_variant_recalibrator->{$mode};
    }

    if ( $self->opts->{'apply_recalibration'} eq 'Y' ) {
        $slurm_jobid_apply_recalibration = CCR::Utilities::General::create_submit_slurm_job(
            debug      => $self->opts->{'debug'},
            account    => $self->opts->{'account'},
            email      => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
            sendemail  => $self->opts->{'sendemail'},
            partition  => $self->opts->{'partition'},
            qos        => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
            cluster    => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
            script     => $scriptname_apply_recalibration,
            command    => $commandline,
            script_dir => join( '/', $self->phaseIII_root, $self->dir->{'phaseIII'}->{'script_apply_recalibration'} ),
            modules    => "java/1.8.0_45,R,gatk/" . $self->opts->{'gatk_version'},
            job_name   => $jobnames_apply_recalibration,
            time       => $gatk_options->{'time'},
            nodes      => 1,
            memory     => $gatk_options->{'mem'},
            ntasks_per_node => $gatk_options->{'nt'},
            dependency      => $dependency,
            logfile          => join( '/',
                $self->phaseIII_root, $self->dir->{'phaseIII'}->{'log_apply_recalibration'},
                $jobnames_apply_recalibration )
        );
    }
    else {
        print "ApplyRecalibration will not be performed!!\n";
    }
    $output_apply_recalibration->{$mode} = $output_file;
    $jobid_apply_recalibration->{$mode}  = $slurm_jobid_apply_recalibration;
    $self->set_recalibrated_vcf($output_apply_recalibration);
    $self->set_jobid_apply_recalibration($jobid_apply_recalibration);
    $self;
}

# merge the bam files created by haplotype caller for exmaining variants under IGV
sub merge_hc_bams {
    my ($self) = shift;
    my $chromosome = $self->get_chromosomes->{'chromosomes_ordered'};

    foreach my $donor ( sort keys %{ $self->config } ) {
        my $sample_id = $self->config->{$donor}->{'case'}->{'ID'};
        my @hc_bams;
        my @dependency;

        foreach my $chr ( @{$chromosome} ) {
            push @hc_bams,    $self->output_hc_bams->{$sample_id}->{$chr};
            push @dependency, $self->jobids_hc->{$sample_id}->{$chr};
        }
        my $input_bam = join( ' ', @hc_bams );
        my $dependency = 'none';

        if ( $self->opts->{'run_haplotype_caller'} eq 'Y' ) {
            $dependency = join( ':', @dependency );
        }

        my $output_dir = join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'hc_bam'} );
        my $output_bam = join( '/', $output_dir,         $self->opts->{'run_id'} . '_' . $donor . '_merged_hc.bam' );

        my $jobname_merge_hc_bams    = join( '_', 'merge_hc_bams', $donor, $self->opts->{'run_id'} );
        my $scriptname_merge_hc_bams = 'slurm_' . $jobname_merge_hc_bams . '.sh';
        my $commandline              = join( ' ', 'samtools merge', $output_bam, $input_bam );
        if ( $self->opts->{'merge_hc_bams'} eq 'Y' ) {
            my $slurm_jobid_merge_hc_bams = CCR::Utilities::General::create_submit_slurm_job(
                debug     => $self->opts->{'debug'},
                email     => $self->opts->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
                sendemail => $self->opts->{'sendemail'},
                partition => $self->opts->{'partition'},
                qos       => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'qos'},
                cluster   => $self->opts->{'partition'} eq 'industry' ? 'industry' : $self->opts->{'cluster'},
                script          => $scriptname_merge_hc_bams,
                command         => $commandline,
                script_dir      => join( '/', $self->phaseII_root, $self->dir->{'phaseII'}->{'script_merge_hc_bams'} ),
                modules         => 'samtools',
                job_name        => $jobname_merge_hc_bams,
                time            => '05:00:00',
                nodes           => 1,
                memory          => 5000,
                account         => $self->opts->{'account'},
                ntasks_per_node => 1,
                dependency      => $dependency,
                logfile          => join( '/',
                    $self->phaseII_root, $self->dir->{'phaseII'}->{'log_merge_hc_bams'},
                    $jobname_merge_hc_bams )
            );
            print "$commandline\n";
        }
    }
}

=head1 AUTHOR

Jianxin Wang, C<< <jw24 at buffalo.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-ngs-tools-gatk at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=NGS-Tools-GATK>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc NGS::Tools::GATK


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=NGS-Tools-GATK>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/NGS-Tools-GATK>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/NGS-Tools-GATK>

=item * Search CPAN

L<http://search.cpan.org/dist/NGS-Tools-GATK/>

=back

=head1 LICENSE AND COPYRIGHT

Copyright 2017 Jianxin Wang.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

no Moose;

__PACKAGE__->meta->make_immutable;

1;    # End of NGS::Tools::GATK

