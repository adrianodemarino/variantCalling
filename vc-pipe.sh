#!/bin/bash

# Arguments

stage=1
slurm=false
job_prefix="${job_prefix}"

while getopts "b:t:f:g:e:p:s:j:" flag; do
case "$flag" in
    b) bed_target=$OPTARG;; ## this is a bed suitable for coverage calculation
    t) target_mutation=$OPTARG;; ## this a bed suitable for target mutation coverage calculation
    f) bed_offtarget=$OPTARG;; ## this a bed suitable for target mutation coverage calculation
    g) simg_dir=$OPTARG;;
    e) stage=$OPTARG;;
    p) padded_bed_target=$OPTARG;; ## this is a bed suitable for variant calling
    s) slurm=$OPTARG;; ## this is a bed suitable for variant calling
    j) job_prefix=$OPTARG;; ## this is a bed suitable for variant calling

esac
done


# Positional arguments
sample_name=${@:$OPTIND:1}
ref=${@:$OPTIND+1:1}
out_dir=${@:$OPTIND+2:1}
pairs1=${@:$OPTIND+3:1}
pairs2=${@:$OPTIND+4:1}


## Auxiliar functions
log(){
     echo "[$(date --rfc-3339=seconds)]: $*"
}

create_cli(){
	sh_cli_fn=$1
	sh_script_fn=$2
	echo $sh_cli_fn > $sh_script_fn
	chmod u+x $sh_script_fn
	job_id_fn=`$sh_script_fn | cut -d" " -f 4`
	echo $job_id_fn
}

# singularity image
if [ "$simg_dir" == "" ]; then
	simg_dir="/data/resources/img/CGT.img"
fi

# get scripts dir
scripts_dir="do_sex_detection_for_future"

# Starting pipeline
log "Sample name:         $sample_name"
log "Reference genome:    $ref"
log "Output directory:    $out_dir"
log "pair1 files:         $pairs1"
log "pair2 files:         $pairs2"
log "bed capture:         $bed_target"
log "Padded bed capture:  $padded_bed_target"
log "Off target bed:      $bed_offtarget"
log "target mutation:     $target_mutation"
log "Singularity dir:     $simg_dir"
log "Scripts dir:         $scripts_dir"
log "Slurm:"	          $slurm

#Example:
#exec_file 		-->			~/repositories/illumina-pipeline/pipeline.sh 

#flag:
#-b 			-->		 	/data/research/projects/veritasgenetics_china/my-bucket-spain/WES_analysis/MedExome_hg19_capture_targets_nochr.bed 
#-t 			-->			/users/jjimenez/repositories/CGT/data/CGT600/v1/variants_CGT600_coverageAnalysis_v1.bed 

#positional argument:
#sample_name 	-->			18-VEX-1009 
#ref 			-->			/data/resources/genome/hg38.fa
#out_dir		--> 		/users/jjimenez/projects/tests_pipelines/analysis/18-VEX-1009 
#pairs1 		--> 		/users/jjimenez/projects/tests_pipelines/rawdata/18-VEX-1009_S2_L001_R1_001.fastq.gz 
#pairs2 		--> 		/users/jjimenez/projects/tests_pipelines/rawdata/18-VEX-1009_S2_L001_R2_001.fastq.gz

# log dir
log_dir=$out_dir/log
mkdir -p $log_dir

# cli dir
cli_dir=$out_dir/cli
mkdir -p $cli_dir

# fastq files
pairs1_array=(`echo $pairs1 | sed 's/,/ /g'`)
pairs2_array=(`echo $pairs2 | sed 's/,/ /g'`)

# checking number of files
if [ ${#pairs1_array[@]} -ne ${#pairs2_array[@]} ]; then
	log "[ERROR] The number of files for pair1 and pair2 must be the same (${#pairs1_array[@]} vs ${#pairs2_array[@]})"
	exit
fi


# Bind path for singularity
bind_path="`dirname $ref`,`dirname $out_dir`"


log "Starting pipeline....................."



# outdir
log "Creating directory $out_dir"
mkdir -p $out_dir


# lanes
lanes=${#pairs1_array[@]}
log "Number of lanes: $lanes"


# Alignment
#log "Starting alignment"

out_bam=$out_dir/${sample_name}

if [ "$stage" -eq 1 ]; then
	log "Doing stage "$stage
	if [ "$lanes" -eq 1 ]; then
		# bam output
		out_bam=$out_dir/${sample_name}

		# bind path of fastq files
		bind_path="$bind_path,`dirname $pairs1`,`dirname $pairs2`"
		log "Singularity bindpath $bind_path"
		job_args=""
		if [ "$slurm" = true ]; then

			sh_cli="sbatch -J ${job_prefix}aln_${sample_name} -o $log_dir/bwa_aln.log --cpus-per-task=6  --mem=10G --wrap=\"singularity exec -B $bind_path $simg_dir bwa mem -M -t 6 -R '@RG\tID:$sample_name\tSM:$sample_name' $ref $pairs1 $pairs2 | singularity exec -B $bind_path $simg_dir samtools view -bSh - -o $out_bam.bam\""
			sh_file="$cli_dir/${job_prefix}aln_${sample_name}.sh"
			job_id=`create_cli "$sh_cli" "$sh_file"`

			sh_cli="sbatch -J ${job_prefix}sort_${sample_name} -o $log_dir/samtools_sort.log --dependency=afterok:$job_id --cpus-per-task=1 --mem=20G --wrap=\"singularity exec -B $bind_path $simg_dir samtools sort $out_bam.bam -o $out_bam.sorted.bam\""
			sh_file="$cli_dir/${job_prefix}sort_${sample_name}.sh"
			job_id=`create_cli "$sh_cli" "$sh_file"`

		else
			singularity exec -B $bind_path $simg_dir bwa mem -M -t 6 -R "@RG\tID:$sample_name\tSM:$sample_name" $ref $pairs1 $pairs2 | singularity exec -B $bind_path $simg_dir samtools view -bSh - | singularity exec -B $bind_path $simg_dir samtools sort -o $out_bam.sorted.bam
		fi


	else
		if [ "$slurm" = true ]; then
			# bam list to merge
			declare -A list_bam
			for index in "${!pairs1_array[@]}"; do
				# lane
				lane=$(($index + 1))
				# fastq files
				pair1=${pairs1_array[$index]}
				pair2=${pairs2_array[$index]}
				# bind path
				bind_path_lane="$bind_path,$(dirname $pair1),$(dirname $pair2)"
				# bam lane
				bam_lane=$out_dir/${sample_name}_lane$lane.bam
				list_bam+="$bam_lane "

				#log "Aligning lane $lane"

				sample_id_lane=${sample_name}_L$lane

				#log "Lane$lane read group id $sample_id_lane"
				#log "Bam output lane1 $bam_lane"

				sh_cli="sbatch -J ${job_prefix}aln${lane}_${sample_name} -o $log_dir/bwa_aln${lane}.log --cpus-per-task=6  --mem=10G --wrap=\"singularity exec -B $bind_path_lane $simg_dir bwa mem -M -t 6 -R '@RG\tID:$sample_id_lane\tSM:$sample_name' $ref $pair1 $pair2 | singularity exec -B $bind_path_lane $simg_dir samtools view -bSh - | singularity exec -B $bind_path_lane $simg_dir samtools sort -o $bam_lane\""
				sh_file="$cli_dir/${job_prefix}aln${lane}_${sample_name}.sh"
				job_id=`create_cli "$sh_cli" "$sh_file"`
				if [ -z "$job_ids" ]; then
				   job_ids="afterok:$job_id"
				else
				   job_ids="$job_ids,afterok:$job_id"
				fi

			done
		fi


		#log "Starting to merge to $out_bam.bam"
		echo $job_ids
		sh_cli="sbatch -J ${job_prefix}merge_${sample_name} -o $log_dir/samtools_merge.log --dependency=$job_ids --cpus-per-task=1 --mem=5G --wrap=\"singularity exec -B $bind_path $simg_dir samtools merge -f $out_bam.sorted.bam ${list_bam[@]} && cp $out_bam.sorted.bam $out_bam.bam && rm ${list_bam[@]} \""
		sh_file="$cli_dir/${job_prefix}merge_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`


	fi
fi # end stage


#log "Starting to sort $out_bam.bam into $out_bam.sorted.bam"
#singularity exec -B $bind_path $simg_dir samtools sort  $out_bam.bam -o $out_bam.sorted.bam
if [ "$stage" -lt 3 ]; then
	log "Doing stage 2"

	#log "Indexing $out_bam.sorted.bam"
	if [ "$slurm" = true ]; then 
		sh_cli="sbatch -J ${job_prefix}idx_${sample_name} -o $log_dir/samtools_idx.log --dependency=afterok:$job_id --cpus-per-task=1 --mem=2G --wrap=\"singularity exec -B $bind_path $simg_dir samtools index $out_bam.sorted.bam\""
		sh_file="$cli_dir/${job_prefix}idx_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`

		#log "Marking duplicates of $out_bam.sorted.bam"
		sh_cli="sbatch -J ${job_prefix}markdup_${sample_name} -o $log_dir/markdup.log --dependency=afterok:$job_id --cpus-per-task=1 --mem=20G --wrap=\"singularity exec -B $bind_path $simg_dir /opt/conda/bin/picard MarkDuplicates INPUT=$out_bam.sorted.bam OUTPUT=$out_bam.sorted.rmDup.bam METRICS_FILE=$out_bam.sorted.rmDup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$out_dir/tmp\""
		sh_file="$cli_dir/${job_prefix}markdup_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`

		sh_cli="sbatch -J ${job_prefix}markdup_index_${sample_name} -o $log_dir/markdup_index.log --dependency=afterok:$job_id --cpus-per-task=1 --mem=20G --wrap=\"singularity exec -B $bind_path $simg_dir samtools index $out_bam.sorted.rmDup.bam\""
		sh_file="$cli_dir/${job_prefix}markdup_index_${sample_name}.sh"
		job_id_rmdupidx=`create_cli "$sh_cli" "$sh_file"`

        #log "Sex detection of $out_bam.sorted.rmDup.bam"
		#sh_cli="sbatch -J ${job_prefix}sex_${sample_name} -o $log_dir/sex.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --mem=20G --wrap=\"$scripts_dir/do_sex_detection.sh $out_bam.sorted.rmDup.bam $simg_dir\""
		#sh_file="$cli_dir/${job_prefix}sex_${sample_name}.sh"
		#job_id=`create_cli "$sh_cli" "$sh_file"`

		#log "Flagstat of $out_bam.sorted.rmDup.bam"
		sh_cli="sbatch -J ${job_prefix}flagstat_${sample_name} -o $log_dir/flagstat.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --mem=20G --wrap=\"singularity exec -B $bind_path $simg_dir samtools flagstat $out_bam.sorted.rmDup.bam > $out_bam.sorted.rmDup.flagstat\""
		sh_file="$cli_dir/${job_prefix}flagstat_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`

		#log "Collect insert size metrics of $out_bam.sorted.rmDup.bam"
		sh_cli="sbatch -J ${job_prefix}isize_${sample_name} -o $log_dir/isize.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --mem=20G --wrap=\"singularity exec -B $bind_path $simg_dir /opt/conda/bin/picard CollectInsertSizeMetrics INPUT=$out_bam.sorted.rmDup.bam OUTPUT=$out_bam.sorted.rmDup.insert.txt H=$out_bam.sorted.rmDup.insert.pdf M=0.5 ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$out_dir/tmp\""
		sh_file="$cli_dir/${job_prefix}isize_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`

		#log "Removig intermediate files"
		sh_cli="sbatch -J ${job_prefix}rm_${sample_name} -o $log_dir/rm.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --wrap=\"rm $out_bam.bam $out_bam.sorted.bam $out_bam.sorted.bam.bai\""
		sh_file="$cli_dir/${job_prefix}rm_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`
	fi

fi

if [ "$stage" -lt 4 ]; then
	log "Doing stage 3"
	genome_file=`echo $ref | sed 's/.fasta$/.genome/'`
	
	if [ "$slurm" = true ]; then
		# create temp bam for coverage
		if [ "$bed_target" != "" ] || [ "$target_mutation" != "" ] || [ "$bed_offtarget" != "" ]; then
			#log "Creating temporal bam for coverage calculation $out_bam.sorted.rmDup.removed.bam"
			sh_cli="sbatch -J ${job_prefix}covbam_${sample_name} -o $log_dir/covbam.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --mem=10G --wrap=\"singularity exec -B $bind_path $simg_dir samtools view -hb -F 1024 -q 1 $out_bam.sorted.rmDup.bam > $out_bam.sorted.rmDup.removed.bam\""
			sh_file="$cli_dir/${job_prefix}covbam_${sample_name}.sh"
			job_id_removed=`create_cli "$sh_cli" "$sh_file"`
		fi

		# calculations of coverage
		if [ "$bed_target" != "" ]; then
        	bind_path="$bind_path,`dirname $bed_target`"
			bed_name=`basename $bed_target .bed`
			out_bed_hist=$out_dir/${sample_name}.${bed_name}.hist.bed
			#log "Histogram coverage calculation of capture bed $bed_target. Output file: $out_bed_hist"
			sh_cli="sbatch -J ${job_prefix}covhist_${sample_name} -o $log_dir/covhist.log --dependency=afterok:$job_id_removed --cpus-per-task=1 --mem=10G --wrap=\"singularity exec -B $bind_path $simg_dir bedtools coverage -sorted -hist -g $genome_file -a $bed_target -b $out_bam.sorted.rmDup.removed.bam | grep ^all > $out_bed_hist\""
			sh_file="$cli_dir/${job_prefix}covhist_${sample_name}.sh"
			job_id_bed_target_hist=`create_cli "$sh_cli" "$sh_file"`
		
			out_bed=$out_dir/${sample_name}.${bed_name}.bed
			#log "Coverage calculation of capture bed $bed_target. Output file: $out_bed"
			sh_cli="sbatch -J ${job_prefix}covbed_${sample_name} -o $log_dir/covbed.log --dependency=afterok:$job_id_removed --cpus-per-task=1 --mem=10G --wrap=\"singularity exec -B $bind_path $simg_dir bedtools coverage -sorted -g $genome_file -a $bed_target -b $out_bam.sorted.rmDup.removed.bam > $out_bed\""
			sh_file="$cli_dir/${job_prefix}covbed_${sample_name}.sh"
			job_id_bed_target=`create_cli "$sh_cli" "$sh_file"`

			dependency_jobs=($job_id_bed_target_hist $job_id_bed_target) # future dependency for removing removed.bam
	
		fi

        # calculation of off target
        if [ "$bed_offtarget" != "" ]; then
            bind_path="$bind_path,`dirname $bed_offtarget`"
	    target_name=`basename $bed_offtarget .bed`
	    out_target=$out_dir/${sample_name}.${target_name}.bed
            #log "Coverage calculation of off target mutation bed $target_mutation. Output file: $out_target"
	    sh_cli="sbatch -J ${job_prefix}cov_offtarget_${sample_name} -o $log_dir/cov_offtarget.log --dependency=afterok:$job_id_removed --cpus-per-task=1 --mem=10G --wrap=\"singularity exec -B $bind_path $simg_dir bedtools coverage -sorted -g $genome_file -a $bed_offtarget -b $out_bam.sorted.rmDup.removed.bam > $out_target\""
	    sh_file="$cli_dir/${job_prefix}cov_offtarget_${sample_name}.sh"
	    job_id_bed_offtarget=`create_cli "$sh_cli" "$sh_file"`
            dependency_jobs+=($job_id_bed_offtarget) # future dependency for removing removed.bam

	fi

	if [ "$target_mutation" != "" ]; then
	    bind_path="$bind_path,`dirname $target_mutation`"
	    target_name=`basename $target_mutation .bed`
	    out_target=$out_dir/${sample_name}.${target_name}.bed
            #log "Coverage calculation of target mutation bed $target_mutation. Output file: $out_target"
	    sh_cli="sbatch -J ${job_prefix}covtarget_${sample_name} -o $log_dir/covtarget.log --dependency=afterok:$job_id_removed --cpus-per-task=1 --mem=10G --wrap=\"singularity exec -B $bind_path $simg_dir bedtools coverage -sorted -g $genome_file -a $target_mutation -b $out_bam.sorted.rmDup.removed.bam > $out_target\""
	    sh_file="$cli_dir/${job_prefix}covtarget_${sample_name}.sh"
	    job_id_target_mutation=`create_cli "$sh_cli" "$sh_file"`

	    dependency_jobs+=($job_id_target_mutation) # future dependency for removing removed.bam
	
	fi

	# remove temp bam created before for coverage calculation
	# create temp bam for coverage
	if [[ -v dependency_jobs ]]; then
            dependency_jobs=$(echo "${dependency_jobs[@]}" | sed 's/ /:/g')
	    #log "Remove $out_bam.sorted.rmDup.removed.bam" 
	    sh_cli="sbatch -J ${job_prefix}covrm_${sample_name} -o $log_dir/covrm.log --dependency=afterok:$dependency_jobs --cpus-per-task=1 --mem=2G --wrap=\"rm $out_bam.sorted.rmDup.removed.bam\""
	    sh_file="$cli_dir/${job_prefix}covrm_${sample_name}.sh"
	    job_id=`create_cli "$sh_cli" "$sh_file"`
	fi
    fi
fi

if [ "$stage" -lt 5 ]; then
	log "Doing stage 4"
	if [ "$slurm" = true ]; then

		if [ "$padded_bed_target" != "" ]; then
			#log "Freebayes of $out_bam.sorted.rmDup.bam with capture bed file $padded_bed_target. Output: $out_bam.vcf"
			sh_cli="sbatch -J ${job_prefix}free_${sample_name} -o $log_dir/free.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --mem=8G --wrap=\"singularity exec -B $bind_path $simg_dir freebayes --min-base-quality 3 -F 0.15 -f $ref -t $padded_bed_target -b $out_bam.sorted.rmDup.bam -v $out_bam.vcf\""
			sh_file="$cli_dir/${job_prefix}free_${sample_name}.sh"
			job_id=`create_cli "$sh_cli" "$sh_file"`
		else
			#log "Freebayes of $out_bam.sorted.rmDup.bam with no capture bed file. Output: $out_bam.vcf"
			sh_cli="sbatch -J ${job_prefix}free_${sample_name} -o $log_dir/free.log --dependency=afterok:$job_id_rmdupidx --cpus-per-task=1 --mem=8G --wrap=\"singularity exec -B $bind_path $simg_dir freebayes --min-base-quality 3 -F 0.15 -f $ref -b $out_bam.sorted.rmDup.bam -v $out_bam.vcf\""
			sh_file="$cli_dir/${job_prefix}free_${sample_name}.sh"
			job_id=`create_cli "$sh_cli" "$sh_file"`
		fi

		sh_cli="sbatch -J ${job_prefix}tabix_${sample_name} -o $log_dir/tabix.log --dependency=afterok:$job_id --cpus-per-task=1 --mem=8G --wrap=\"singularity exec -B $bind_path $simg_dir bgzip $out_bam.vcf && singularity exec -B $bind_path $simg_dir tabix -p vcf $out_bam.vcf.gz\""
		sh_file="$cli_dir/${job_prefix}bgzip_${sample_name}.sh"
		job_id=`create_cli "$sh_cli" "$sh_file"`
		#log "$job_id"
	fi
fi
