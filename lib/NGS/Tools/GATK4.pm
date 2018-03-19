package NGS::Tools::GATK4;

use strict;
use warnings;
use Moose;
use v5.10;
use Carp;
use File::Basename;
use Params::Validate qw(:all);

use CCR::Utilities::General;

our $VERSION = '1.0.0';
our $jobid_vqsr;
our $output_vqsr;
our $tranches_vqsr;
our $jobid_avqsr;
our $output_avqsr;

has 'bam' => (
	is => 'ro',
	isa => 'HashRef'
);

has 'config' => (
	is => 'ro',
	isa => 'HashRef'
);

has 'opts' => (
	is => 'ro',
	isa => 'HashRef'
);

has 'subdirs' => (
	is => 'rw',
	isa => 'HashRef'
);

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

has 'chromosome' => (
	is  => 'rw',
	isa => 'ArrayRef'
);

has 'intervals' => (
	is  => 'rw',
	isa => 'HashRef'
);
has 'jobid_base_recalibration' => (
    is     => 'rw',
	isa    => 'HashRef'
);

has 'recalibration_table' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'recalibrated_bams' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'jobids_print_recalibrated_bam' => (
	  is => 'rw',
    	isa => 'HashRef'
);

has 'gvcfs' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'jobids_hc' => (
	is => 'rw',
	isa => 'HashRef'
);


has 'jobids_importGDB' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'gemomicsDBs' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'jobids_genotypeGVCF' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'raw_vcfs' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'jobid_gather_vcfs' => (
	is => 'rw',
	isa => 'Str'
);

has 'output_gather_vcfs' => (
	is => 'rw',
	isa => 'Str'
);

has 'jobid_recal' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'output_recal' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'output_tranches' => (
	is => 'rw',
	isa => 'HashRef'
);

has 'jobid_apply_vqsr' => (
	 is => 'rw',
     isa => 'HashRef'
);

has 'output_apply_vqsr' => (
      is => 'rw',
      isa => 'HashRef'
);


# define subdirectory structure
sub dir_structure {
    my $self = shift;

    my %dirs = (
        $self->phaseI_root => {
            intermediate_bam_files         => 'intermediate_bam_files',
			intermed_batch_id              => 'intermediate_bam_files/' .$self->opts->{'batch_id'},
            recalibration_table            => 'recalibration_table',
            scripts                        => 'scripts',
            script_recalibration_table     => 'scripts/recalibration_table',
            script_create_recalibrated_bam => 'scripts/create_recalibrated_bam',
			script_cr_recal_bam_batchid      => 'scripts/create_recalibrated_bam/' .$self->opts->{'batch_id'},
            log                            => 'log',
            log_recalibration_table        => 'log/recalibration_table',
            log_create_recalibrated_bam    => 'log/create_recalibrated_bam',
			log_cr_recal_bam_batchid         => 'log/create_recalibrated_bam/' .$self->opts->{'batch_id'},
        },
        $self->phaseII_root => {
            raw_vcf                        => 'raw_vcf',
            scripts                        => 'scripts',
            script_joint_calling_genotypes => 'scripts/joint_calling_genotypes',
			script_joint_cal_gt_batchid      => 'scripts/joint_calling_genotypes/' .$self->opts->{'batch_id'},
            script_haplotype_caller        => 'scripts/haplotype_caller',
			script_haplotype_caller_batchid  => 'scripts/haplotype_caller/' .$self->opts->{'batch_id'},
			script_importgenomicsdb   => 'scripts/importGenimcsDB',
            script_gather_vcfs            => 'scripts/gather_vcfs',
            script_merge_hc_bams           => 'scripts/merge_hc_bams',
     		script_merge_hc_bams_batchid     => 'scripts/merge_hc_bams/' .$self->opts->{'batch_id'},
            hc_bam                         => 'hc_bam',
			gemomicsDB       => 'genomicsDB',
            log                            => 'log',
			log_importgenomicsdb       => 'log/importGenomicsDB',
            log_joint_calling_genotypes    => 'log/joint_calling_genotypes',
			log_joint_calling_gt_batchid     => 'log/joint_calling_genotypes/' .$self->opts->{'batch_id'},
            log_haplotype_caller           => 'log/haplotype_caller',
			log_haplotype_caller_batchid    => 'log/haplotype_caller/' .$self->opts->{'batch_id'},
            log_gather_vcfs               => 'log/gather_vcfs',
            log_merge_hc_bams              => 'log/merge_hc_bams',
			log_merge_hc_bams_batchid        => 'log/merge_hc_bams/' .$self->opts->{'batch_id'},
        },
        $self->phaseIII_root => {
            scripts                     => 'scripts',
            recalibrated_vcf            => 'recalibrated_vcf',
            vcf_recalibration_model     => 'vcf_recalibration_model',
            log                         => 'log',
        }
    );

    $self->subdirs( \%dirs );
    $self;
}

sub run_base_recalibration {
    my $self = shift;
    my $jobids_base_recalib;
    my $recal_tables;

    foreach my $donor ( sort keys %{ $self->bam } ) {
		foreach my $sample_type (qw(case control)) {
	        my $sample_id   = $self->bam->{$donor}->{$sample_type}->{'ID'};
			next if (! $sample_id); # the specific $smaple_type does existi (i.e. not a case/control or tumor/normal setup), so skip
	        my $input_bam   = $self->bam->{$donor}->{$sample_type}->{'path'};
	        my $output_file = join( '_', $self->opts->{'batch_id'}, $sample_id, 'recalibration_table.grp' );
	        my $output_dir  = join( '/', $self->phaseI_root, $self->subdirs->{$self->phaseI_root}->{recalibration_table});
	        $output_file = join( '/', $output_dir, $output_file );
	
				
     		my $known_sites = '--known-sites '. join( ' --known-sites ',
            $self->config->{'dbsnp'},  $self->config->{'g1k_snp'}, $self->config->{'indels'},
            $self->config->{'hapmap'}, $self->config->{'mills'} );

			# java need to have Xmx values slightly less than the actual memory required by slurm
			my $xmx = $self->config->{BaseRecalibrator}->{mem} - 2000; # make it 2gb smaller
   			my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'BaseRecalibrator', '-I', $input_bam,
							    '-R', $self->config->{ref}, $known_sites, '-O', $output_file,
								'--verbosity', 'INFO',
 								'--TMP_DIR', $self->config->{TMP_DIR});
	
	        my $jobname_base_recalibration = join( '_', 'BQR', $sample_id, $self->opts->{'batch_id'} );
	        my $scriptname_base_recalibration = 'slurm_' . $jobname_base_recalibration . ".sh";
	        my $slurm_jobid_crt = 'none';

	        if ( $self->config->{'BaseRecalibrator'}->{'run'} eq 'Yes' ) {
	            $slurm_jobid_crt = CCR::Utilities::General::create_submit_slurm_job(
	                debug      => $self->config->{'debug'},
	                account    => $self->config->{'account'},
	                email      => $self->config->{'email'} || $self->opts->{'user'}. '@buffalo.edu',
	                sendemail  => $self->config->{'sendemail'},
	                partition  => $self->config->{'partition'},
	                qos        => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
	                cluster    => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
	                script     => $scriptname_base_recalibration,
	                command    => $command,
	                script_dir => $self->phaseI_root .'/'. $self->subdirs->{$self->phaseI_root}->{'script_recalibration_table'},
	                modules    => "java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
	                job_name   => $jobname_base_recalibration,
	                time       => $self->config->{BaseRecalibrator}->{'time'},
	                nodes      => 1,
	                memory     => $self->config->{BaseRecalibrator}->{'mem'},
	                ntasks_per_node => $self->config->{BaseRecalibrator}->{'nct'},
	                dependency      => 'none',
	                logfile          => $self->phaseI_root .'/'. $self->subdirs->{$self->phaseI_root}->{'log_recalibration_table'} .'/'. $jobname_base_recalibration
	            );
	
	            $jobids_base_recalib->{$sample_id} = $slurm_jobid_crt;
	        }
	        else {
	            print "Skipping BaseRecalibrator...\n";
	        }
	
	        $recal_tables->{$sample_id} = $output_file;
	    	$self->jobid_base_recalibration($jobids_base_recalib) if ($jobids_base_recalib);
	    	$self->recalibration_table($recal_tables);
		}
	}
	$self;
}

sub create_recalibrated_bams{
	my $self = shift;

    my $jobids_create_recalibrated_bam;
    my $jobname_create_recalibrated_bam;
    my $recalibrated_bams;

    foreach my $donor ( sort keys %{ $self->bam } ) {
		foreach my $sample_type (qw(case control)) {
	        my $sample_id   = $self->bam->{$donor}->{$sample_type}->{'ID'};
			next if (! $sample_id); # the specific $smaple_type does existi (i.e. not a case/control or tumor/normal setup), so skip
	        my $input_bam   = $self->bam->{$donor}->{$sample_type}->{'path'};

	        my $output_dir  = join( '/', $self->phaseI_root, $self->subdirs->{$self->phaseI_root}->{intermed_batch_id});
	        my $output_file = join( '_', $self->opts->{'batch_id'}, $sample_id, 'recalibrated.bam' );
	        $output_file = join( '/', $output_dir, $output_file );

			# java need to have Xmx values slightly less than the actual memory required by slurm
            my $xmx = $self->config->{ApplyBQSR}->{mem} - 2000; # make it 2gb smaller
            my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'ApplyBQSR', '-I', $input_bam,
                                '-R', $self->config->{ref}, '--bqsr-recal-file', $self->recalibration_table->{$sample_id},
								'-O', $output_file,
                                '--verbosity', 'INFO',
                                '--TMP_DIR', $self->config->{TMP_DIR});

			my $scriptname_create_recalibrated_bam =
		  		join( '_', 'slurm', 'BQRB', $sample_id, $self->opts->{'batch_id'} . "\.sh" );
			my $jobname_create_recalibrated_bam = join( '_', 'BQRB', $sample_id, $self->opts->{'batch_id'} );

			my $slurm_jobid_create_recalibrated_bam = 'none';
			my $dependency = 'none';
			if ( $self->config->{'ApplyBQSR'}->{'run'} eq 'Yes' ) {
				$dependency = $self->jobid_base_recalibration->{$sample_id};
			}

			if ( $self->config->{'ApplyBQSR'}->{'run'} eq 'Yes' ) {
				$slurm_jobid_create_recalibrated_bam = CCR::Utilities::General::create_submit_slurm_job(
					debug     => $self->config->{'debug'},
					account   => $self->config->{'account'},
					email     => $self->config->{'email'} || $self->config->{'user'}. '@buffalo.edu',
					sendemail => $self->config->{'sendemail'},
					partition => $self->config->{'partition'},
					qos       => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
					cluster   => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
					script    => $scriptname_create_recalibrated_bam,
					command   => $command,
					script_dir =>
				  		join( '/', $self->phaseI_root, $self->subdirs->{$self->phaseI_root}->{'script_cr_recal_bam_batchid'} ),
					modules         => "java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
					job_name        => $jobname_create_recalibrated_bam,
					time            => $self->config->{'ApplyBQSR'}->{'time'},
					nodes           => 1,
					memory          => $self->config->{'ApplyBQSR'}->{'mem'},                                  #12000,
					ntasks_per_node => $self->config->{'ApplyBQSR'}->{'nct'},                                  #4,
					dependency      => $dependency,
					logfile          => join( '/',
						$self->phaseI_root, $self->subdirs->{$self->phaseI_root}->{'log_cr_recal_bam_batchid'},
						$jobname_create_recalibrated_bam )
				);
			}
			else {
				print "Skipping BQSR...\n";
			}

			$recalibrated_bams->{$sample_id} = $output_file;

			$jobids_create_recalibrated_bam->{$sample_id} = $slurm_jobid_create_recalibrated_bam
		  	if ( defined $slurm_jobid_create_recalibrated_bam );

			$self->recalibrated_bams($recalibrated_bams);
			$self->jobids_print_recalibrated_bam($jobids_create_recalibrated_bam);
		}
	}
    $self;
}

sub run_haplotype_caller {
    my $self = shift;
    my $jobids_haplotype_caller;
    my $output_gvcfs;

    foreach my $chr ( @{$self->chromosome} ) {
		my $interval_list = $self->config->{interval_list};
		$interval_list =~ s/regions/regions_$chr/;

		# create a dir to store the output file
		my $output_dir = join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'raw_vcf'}, $chr );
		File::Path::make_path( File::Spec->rel2abs($output_dir) );

		my $output_dir_hc_bam = join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'hc_bam'}, $chr);
		File::Path::make_path( File::Spec->rel2abs($output_dir_hc_bam) );

		foreach my $donor ( sort keys %{ $self->bam } ) {
        	foreach my $sample_type (qw(case control)) {
            	my $sample_id   = $self->bam->{$donor}->{$sample_type}->{'ID'};
            	next if (! $sample_id); # the specific $smaple_type does existi (i.e. not a case/control or tumor/normal setup), so skip
            	my $input_bam   = $self->bam->{$donor}->{$sample_type}->{'path'};
	            my $output_file = join( '_', $sample_id, $chr, 'g.vcf.gz' );
    	        $output_file = join( '/', $output_dir, $output_file );
			
				my $hc_realigned_bam = join( '_', $sample_id, $chr, 'HC.bam');
				$hc_realigned_bam = join( '/', $output_dir_hc_bam,  $hc_realigned_bam);

            	# java need to have Xmx values slightly less than the actual memory required by slurm
            	my $xmx = $self->config->{HaplotypeCaller}->{mem} - 2000; # make it 2gb smaller
            	my $command = join( ' ', 'gatk --java-options', "\"-Xmx${xmx}M\"", 'HaplotypeCaller -I', $input_bam,
                                '-R', $self->config->{ref}, '-L', $interval_list,
                                '-O', $output_file, '-bamout', $hc_realigned_bam, '-ERC GVCF',
                                '--verbosity INFO',
                                '--TMP_DIR', $self->config->{TMP_DIR});

				if ($self->config->{'DataType'} eq 'RNA-Seq'){
					$command .= ' -dontUseSoftClippedBases -stand_call_conf 20.0';
				}

            	my $jobname_haplotype_caller = join( '_', 'HC', $sample_id, $chr, $self->opts->{'batch_id'} );
				my $scriptname_haplotype_caller = 'slurm_' . $jobname_haplotype_caller . "\.sh";

				my $dependency = $self->jobids_print_recalibrated_bam->{$sample_id} || 'none';

				my $slurm_jobid_hc = 'none';
				if ( $self->config->{'HaplotypeCaller'}->{'run'} eq 'Yes' ) {
					$slurm_jobid_hc = CCR::Utilities::General::create_submit_slurm_job(
						debug           => $self->config->{'debug'},
						account         => $self->config->{'account'},
						email           => $self->config->{'email'}|| $self->config->{'user'}. '@buffalo.edu',
						sendemail       => $self->config->{'sendemail'},
						partition       => $self->config->{'partition'},
						qos             => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
						cluster         => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
						script          => $scriptname_haplotype_caller,
						command         => $command,
						script_dir =>
					  			join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'script_haplotype_caller_batchid'} ),
						modules         => "java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
						job_name        => $jobname_haplotype_caller,
						time            => $self->config->{'HaplotypeCaller'}->{'time'},
						nodes           => 1,
						memory          => $self->config->{'HaplotypeCaller'}->{'mem'},
						ntasks_per_node => $self->config->{'HaplotypeCaller'}->{'nct'},
						dependency      => $dependency,
						logfile          => join( '/',
							$self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'log_haplotype_caller_batchid'},
							$jobname_haplotype_caller )
					);
				} else {
					print "Skipping HaplotypeCaller...\n";
				}
				push @{$jobids_haplotype_caller->{$chr}}, $slurm_jobid_hc;
				push @{$output_gvcfs->{$chr}}, $output_file;
			}
	
			$self->gvcfs($output_gvcfs);
			$self->jobids_hc($jobids_haplotype_caller);
        	}
	}
    $self;
}

sub run_GenomicsDBimport {
	my $self = shift;	
	my $jobids_importGDB;
	my $gemomicsdbs;

	my $output_dir = join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'gemomicsDB'} );
	foreach my $chr ( @{$self->chromosome} ) {
		# make a sampel ID to gvcf file path map
		my $map_file = $output_dir .'/'. $chr . '_map.txt';
		open(MP, ">$map_file") or croak "Error open $map_file: $!\n";
		
		foreach my $gvcf (@{$self->gvcfs->{$chr}}) {

			my @tmp = split(/\//, $gvcf);
			
			my ($sample_id) = $tmp[-1] =~ /^(.*?)_chr.*/;
		#	open(GV, "gunzip -c $gvcf |") or croak "Error open $gvcf: $!\n";

		#	while(<GV>){
		#		if ($_ =~ /^\#CHROM/){
		#			chomp;
		#			my @tmp = split(/\t/);
		#			print MP "$tmp[9]\t$gvcf\n";
		#			last;
		#		}
		#	}
		#	close GV;
			print MP "$sample_id\t$gvcf\n";
		}
		close MP;
		
		# loop through each interval and submit slurm jobs
		foreach my $interval (@{$self->intervals->{$chr}}) {
				my $output_file = $output_dir .'/'. 'GenomicDB_'. $self->opts->{'batch_id'} .'_'. $interval;

				my $xmx = $self->config->{GenomicsDBImport}->{mem} - 4000;

				my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'GenomicsDBImport', '--sample-name-map', $map_file, 
					'--genomicsdb-workspace-path', $output_file, '--overwrite-existing-genomicsdb-workspace true',
					'-ip 0 --batch-size 50 --reader-threads 5', '-L', $interval);

				my $dependency = join(':', @{$self->jobids_hc->{$chr}}) || 'none';
				$dependency = 'none' if ($dependency =~ /none/);
				my $jobname_importGDB = join('_', 'IGDB', $interval, $self->opts->{batch_id});
				my $script_name_importGDB = 'slurm_' . $jobname_importGDB . '.sh';

				my $slurm_jobid_importGDB = 'none';
				if ($self->config->{'GenomicsDBImport'}->{'run'} eq 'Yes' ){
					$slurm_jobid_importGDB = CCR::Utilities::General::create_submit_slurm_job(
							debug           => $self->config->{'debug'},
							account         => $self->config->{'account'},
							email           => $self->config->{'email'}|| $self->config->{'user'}. '@buffalo.edu',
							sendemail       => $self->config->{'sendemail'},
							partition       => $self->config->{'partition'},
							qos             => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
							cluster         => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
							script          => $script_name_importGDB,
							command         => $command,
							script_dir =>
							  join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'script_importgenomicsdb'} ),
							modules         => "java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
							job_name        => $jobname_importGDB,
							time            => $self->config->{'GenomicsDBImport'}->{'time'},
							nodes           => 1,
							memory          => $self->config->{'GenomicsDBImport'}->{'mem'},
							ntasks_per_node => $self->config->{'GenomicsDBImport'}->{'nct'},
							dependency      => $dependency,
							logfile          => join( '/',
								$self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'log_importgenomicsdb'},
								$jobname_importGDB )
						);
				} else {
					print "Skipping GenomicsDBImport...\n";
				}
				$jobids_importGDB->{$interval} = $slurm_jobid_importGDB;
				$gemomicsdbs->{$interval} = $output_file;
			}
	}
	$self->jobids_importGDB($jobids_importGDB);
	$self->gemomicsDBs($gemomicsdbs);
	$self;
}

sub joint_calling_genotypes {
	my $self = shift;
	my $jobids_gt;
	my $raw_vcfs;
	my $output_dir = join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'raw_vcf'} );	
	
	foreach my $chr (@{$self->chromosome} ) {
		foreach my $interval (@{$self->intervals->{$chr}}) {
			my $gendb = $self->gemomicsDBs->{$interval};
			my $output_vcf = $output_dir .'/'. $interval . '_'. $self->opts->{'batch_id'}. '_raw.vcf.gz';
			my $xmx = $self->config->{GenotypeGVCFs}->{mem} - 2000;
			my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'GenotypeGVCFs', 
			'-R', $self->config->{ref}, '-V', 'gendb://'. $gendb, '-O', $output_vcf, '--dbsnp', $self->config->{dbsnp}, '-L', $interval, '-G StandardAnnotation' );
	
			my $slurm_jobid_jcg = 'none';
			my $dependency = $self->jobids_importGDB->{$interval} || 'none';
			my $jobname_jcg = 'JCG_' . $interval.'_'. $self->opts->{batch_id}; 
			my $script_name_jcg = 'slurm_' . $jobname_jcg . '.sh';

			if ($self->config->{'GenotypeGVCFs'}->{'run'} eq 'Yes' ) {
				$slurm_jobid_jcg = CCR::Utilities::General::create_submit_slurm_job(
					debug           => $self->config->{'debug'},
					account         => $self->config->{'account'},
					email           => $self->config->{'email'}|| $self->config->{'user'}. '@'. $self->config->{emailDomain},
					sendemail       => $self->config->{'sendemail'},
					partition       => $self->config->{'partition'},
					qos             => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
					cluster         => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
					script          => $script_name_jcg,
					command         => $command,
					script_dir =>
					  join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'script_joint_cal_gt_batchid'} ),
					modules         => "java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
					job_name        => $jobname_jcg,
					time            => $self->config->{'GenotypeGVCFs'}->{'time'},
					nodes           => 1,
					memory          => $self->config->{'GenotypeGVCFs'}->{'mem'},
					ntasks_per_node => $self->config->{'GenotypeGVCFs'}->{'nct'},
					dependency      => $dependency,
					logfile          => join( '/',
						$self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'log_joint_calling_gt_batchid'},
						$jobname_jcg )
				);
			} else {
				print "Skipping GenotypeGVCFs...\n";
			}
			$jobids_gt->{$interval} = $slurm_jobid_jcg;
			$raw_vcfs->{$interval} = $output_vcf;
		}
	}
	$self->jobids_genotypeGVCF($jobids_gt);	
	$self->raw_vcfs($raw_vcfs);
	$self;
}

sub gather_vcfs {
    my $self         = shift;
    my @input_vcfs;
    my @dependency;
    foreach my $chr ( @{$self->chromosome} ) {
		foreach my $interval (@{$self->intervals->{$chr}}){
        	push @input_vcfs, $self->raw_vcfs->{$interval};
        	push @dependency,    $self->jobids_genotypeGVCF->{$interval};
    	}
	}
    my $input_vcfs   = join( " -I ", @input_vcfs );
    my $output_dir  = join( '/',    $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'raw_vcf'} );
    my $output_file = join( '/',    $output_dir,  $self->opts->{'batch_id'} .'_merged_raw.vcf.gz' );
	my $xmx = $self->config->{GatherVcfs}->{mem} - 2000; 
	my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'GatherVcfs', '-I', $input_vcfs, '-O', $output_file);
		
    my $jobname_gather_vcfs     = 'GaVCF_'. $self->opts->{'batch_id'};
    my $scriptname_gather_vcfs  = 'slurm_' . $jobname_gather_vcfs .'_'. $self->opts->{'batch_id'} .'.sh';
    my $dependency               = join( ':', @dependency );
	$dependency = 'none' if ($dependency =~ /none/);
    my $slurm_jobid_gather_vcf = 'none';

    if ( $self->config->{'GatherVcfs'}->{'run'} eq 'Yes' ) {
        $slurm_jobid_gather_vcf = CCR::Utilities::General::create_submit_slurm_job(
            debug           => $self->config->{'debug'},
            account         => $self->config->{'account'},
            email           => $self->config->{'email'} || $self->config->{'user'}. '@'. $self->config->{emailDomain},
            sendemail       => $self->config->{'sendemail'},
            partition       => $self->config->{'partition'},
            qos             => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
            cluster         => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
            script          => $scriptname_gather_vcfs,
            command         => $command,
            script_dir      => join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'script_gather_vcfs'} ),
            modules         => "java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
            job_name        => $jobname_gather_vcfs,
            time            => $self->config->{'GatherVcfs'}->{'time'},
            nodes           => 1,
            memory          => $self->config->{'GatherVcfs'}->{'mem'},
            ntasks_per_node => $self->config->{'GatherVcfs'}->{'nt'},
            dependency      => $dependency,
            logfile         =>
              join( '/', $self->phaseII_root, $self->subdirs->{$self->phaseII_root}->{'log_gather_vcfs'}, $jobname_gather_vcfs )
        );
    } else {
		print "Skipping GatherVcfs...\n";
	}
    $self->jobid_gather_vcfs($slurm_jobid_gather_vcf);
    $self->output_gather_vcfs($output_file);
    $self;
}

sub run_variant_recalibrator {
	my $self = shift;
	my $mode = shift;

 	# ensure the correct values are passed
 	if ( $mode !~ /SNP|INDEL/){
		croak "Invalide values passed. Valide value is 'SNP' or 'INDEL'";
	}

	my $output_dir  = join( '/',    $self->phaseIII_root, $self->subdirs->{$self->phaseIII_root}->{'vcf_recalibration_model'} );
	my $output_recal = join( '/',    $output_dir, $self->opts->{'batch_id'} . '_recalibration_'. $mode );
	my $output_tranches = join( '/',    $output_dir, $self->opts->{'batch_id'} . '_tranches_'. $mode );
	my $output_rscript = join( '/',    $output_dir, $self->opts->{'batch_id'} . '.Rscript');
	my $xmx = $self->config->{'VariantRecalibrator'}->{'mem'} - 1000;
	my $resource;

	if ($mode eq 'SNP') {
		$resource =  join(" ", 
   						'--resource hapmap,known=false,training=true,truth=true,prior=15.0:' . $self->config->{hapmap},
   						'--resource omni,known=false,training=true,truth=false,prior=12.0:'. $self->config->{omni},
   						'--resource 1000G,known=false,training=true,truth=false,prior=10.0:'. $self->config->{g1k_snp},
   						'--resource dbsnp,known=true,training=false,truth=false,prior=2.0:'. $self->config->{dbsnp},
						);
	} elsif ($mode eq 'INDEL'){
		$resource =  join(" ",
					 '--resource mills,known=false,training=true,truth=true,prior=12.0:'. $self->config->{'mills'},
                	'--resource dbsnp,known=true,training=false,truth=false,prior=2.0:'. $self->config->{'dbsnp'},
 					'--max-gaussians 4');

	}
	my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'VariantRecalibrator', 
						'-R', $self->config->{'ref'}, '-V', $self->output_gather_vcfs, $resource,
						'-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR',
						'-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.5 -tranche 99.3 -tranche 99.0 -tranche 90.0',
   						'--output',  $output_recal,
   						'--tranches-file', $output_tranches,
   						'--rscript-file',  $output_rscript
						);
	$command =~ s/\-an QD// if ($self->config->{DataType} eq 'WES');
	my $jobname_recal = 'VRecal_'. $mode .'_'. $self->opts->{'batch_id'};
	my $scriptname_recal = 'slurm_'. $jobname_recal .'.sh';
	my $dependency = $self->jobid_gather_vcfs || 'none';
	my $slurm_jobid_recal = 'none';

    if ( $self->config->{'VariantRecalibrator'}->{'run'} eq 'Yes' ) {
        $slurm_jobid_recal = CCR::Utilities::General::create_submit_slurm_job(
            debug           => $self->config->{'debug'},
            account         => $self->config->{'account'},
            email           => $self->config->{'email'} || $self->config->{'user'}. '@'. $self->config->{emailDomain},
            sendemail       => $self->config->{'sendemail'},
            partition       => $self->config->{'partition'},
            qos             => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
            cluster         => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
            script          => $scriptname_recal,
            command         => $command,
            script_dir      => join( '/', $self->phaseIII_root, $self->subdirs->{$self->phaseIII_root}->{'scripts'} ),
            modules         => "R,java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
            job_name        => $jobname_recal,
            time            => $self->config->{'VariantRecalibrator'}->{'time'},
            nodes           => 1,
            memory          => $self->config->{'VariantRecalibrator'}->{'mem'},
            ntasks_per_node => $self->config->{'VariantRecalibrator'}->{'nt'},
            dependency      => $dependency,
            logfile         =>
              join( '/', $self->phaseIII_root, $self->subdirs->{$self->phaseIII_root}->{'log'}, $jobname_recal )
        );
    } else {
		print "Skipping VariantRecalibrator...\n";
	}
	$jobid_vqsr->{$mode} = $slurm_jobid_recal;
	$output_vqsr->{$mode} = $output_recal;
	$tranches_vqsr->{$mode} = $output_tranches;

    $self->jobid_recal($jobid_vqsr);
    $self->output_recal($output_vqsr);
	$self->output_tranches($tranches_vqsr);
    $self;

}

sub run_apply_recalibration {
	my $self = shift;
	my $mode = shift;
	my $output_dir  = join( '/',    $self->phaseIII_root, $self->subdirs->{$self->phaseIII_root}->{'recalibrated_vcf'} );
	my ($input_vcf, $output_vcf);

	if ( $mode eq 'INDEL'){
		$input_vcf = $self->output_apply_vqsr->{'SNP'};
		$output_vcf = $self->output_apply_vqsr->{'SNP'};
		$output_vcf =~ s/recalibrated/INDEL_recalibrated/;
	} else {
		$input_vcf = $self->output_gather_vcfs;
		$output_vcf = join( '/',    $output_dir,  $self->opts->{'batch_id'} .'_apply_recal_'. $mode . '_recalibrated.vcf.gz' );
	}

	my $xmx = $self->config->{ApplyVQSR}->{mem} - 1000;
	my $command = join( ' ', 'gatk', '--java-options', "\"-Xmx${xmx}M\"", 'ApplyVQSR',
						'-R', $self->config->{ref}, '-V', $input_vcf,
						'-O', $output_vcf,
						'--tranches-file', $self->output_tranches->{$mode},
						'--recal-file', $self->output_recal->{$mode}, '-mode', $mode);
	my $dependency = 'none';
	if ($mode eq 'INDEL'){
		$dependency = $self->jobid_apply_vqsr->{'SNP'};
		$command .= ' -ts-filter-level 95.0';
	} else {
		$dependency = $self->jobid_recal->{$mode} || 'none';
		$command .= ' -ts-filter-level 99.0';
	}

	my $jobname_apvqsr = 'AVQSR_' . $mode .'_'. $self->opts->{'batch_id'};
	my $script_name_apvqsr = 'slurm_' . $jobname_apvqsr . '.sh';
	my $slurm_jobid_apvqsr = 'none';

	if ($self->config->{'ApplyVQSR'}->{'run'} eq 'Yes'){	
		$slurm_jobid_apvqsr = CCR::Utilities::General::create_submit_slurm_job(
            debug           => $self->config->{'debug'},
            account         => $self->config->{'account'},
            email           => $self->config->{'email'} || $self->config->{'user'}. '@'. $self->config->{emailDomain},
            sendemail       => $self->config->{'sendemail'},
            partition       => $self->config->{'partition'},
            qos             => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'qos'},
            cluster         => $self->config->{'partition'} eq 'industry' ? 'industry' : $self->config->{'cluster'},
            script          => $script_name_apvqsr,
            command         => $command,
            script_dir      => join( '/', $self->phaseIII_root, $self->subdirs->{$self->phaseIII_root}->{'scripts'} ),
            modules         => "R,java/1.8.0_45,gatk/" . $self->config->{'GATK_version'},
            job_name        => $jobname_apvqsr,
            time            => $self->config->{'VariantRecalibrator'}->{'time'},
            nodes           => 1,
            memory          => $self->config->{'VariantRecalibrator'}->{'mem'},
            ntasks_per_node => $self->config->{'VariantRecalibrator'}->{'nt'},
            dependency      => $dependency,
            logfile         =>
              join( '/', $self->phaseIII_root, $self->subdirs->{$self->phaseIII_root}->{'log'}, $jobname_apvqsr )
        );
	}	
	$jobid_avqsr->{$mode} = $slurm_jobid_apvqsr;
	$output_avqsr->{$mode} = $output_vcf;
	$self->jobid_apply_vqsr($jobid_avqsr);
	$self->output_apply_vqsr($output_avqsr);
	$self;
}

# get chromosome names from input bam files
sub get_chromosomes {
    my $self = shift;

	my $input_bam;	
	foreach my $donor ( sort keys %{ $self->bam } ) {
        foreach my $sample_type (qw(case control)) {
            my $sample_id   = $self->bam->{$donor}->{$sample_type}->{'ID'};
            next if (! $sample_id); # the specific $smaple_type does existi (i.e. not a case/control or tumor/normal setup), so skip
            $input_bam   = $self->bam->{$donor}->{$sample_type}->{'path'};

			# check if file exists. If so, exit the loop (assume all the bams are aligned using the same reference)
			last if (-e $input_bam);
        }
    }

    # get the header section for just one bam file (assuming all bam files are aligned on the same ref genome)
    my $bam_header = `samtools view -H $input_bam`;
	my @chromosomes;

    # extract chromosome/contig names from bam header
    while ($bam_header =~ /\@SQ\s+SN:(.*?)\s+LN:(\d+)/g ) {
		my $chr = $1;
        # skip non chromosome contigs
        next if ($chr =~ /_/);
		next if ($chr !~ /[0-9XY]$/);
     	push @chromosomes, $chr;
    }
	$self->chromosome(\@chromosomes);
}

sub get_intervals {
	my $self = shift;
	my $interval_list = $self->config->{interval_list};
	my $intervals_;
	my $interval_size = 2_000_000; # set to 2MB
	my $tmp;
	
	open(IN, "$interval_list") or croak "Error open $interval_list: $!\n";

	while(<IN>){
		next if (/^\@/);
		my ($chr, $start, $end) = split(/\t/);
		my $range = $end - $start + 1;

		while(1){
			if ($range <= 1.5*$interval_size){
				push @{$intervals_->{$chr}}, $chr .':'. $start .'-'.  $end;
				last;
			
			} else {
				$tmp = $start + $interval_size;
				push @{$intervals_->{$chr}}, $chr .':'. $start .'-'.  $tmp;
				
				$range = $end - $tmp;
				$start = $tmp + 1;
			}
		}
	}
	close IN;
	
	$self->intervals($intervals_);
	$self;

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

