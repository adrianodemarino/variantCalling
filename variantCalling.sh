#!/bin/bash
#pipe for varinat calling
genome=/opt/NGS/genomes/hg38.p12.fa
tmpdir=/opt/NGS/scratch/
biotools=/home/adriano/src/img/biotools.img
mount_dir_singu=/home/adriano/
inputdir=/opt/NGS/data/memorial_hospital
outdir=/opt/NGS/results/memorial_hospital
sample_list=/opt/NGS/data/sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt
for idsample in $(cat $sample_list); 
	
	do 
	idsample_R1=`ls *$idsample*_1*`
	idsample_R2=`ls *$idsample*_2*`
	echo Start analysis of sample $idsample_R1 $idsample_R2
	#1. Dowload human hg38 genome
	#2. index it with: `bwa index -a bwtsw hg38.fa.gz`
	#3. mapping
	echo " mapping bwa "
	bwa mem -t 8 -R "@RG\tID:$idsample\tSM:$idsample"  $inputdir/$idsample_R1 $inputdir/$idsample_R2 | samtools view -b - > $outdir/$idsample.raw.bam
	#4. View
	#singularity exec -B $mount_dir_singu $biotools samtools view -h $outdir/$idsample.raw.bam | less -S
	#5. 
	echo " sambamba sort "
	singularity exec -B $mount_dir_singu $biotools sambamba sort -t 10 -m 25G --tmpdir $tmpdir -o $outdir/${idsample}.raw.sorted.bam $outdir/${idsample}.raw.bam
	#6. remove duplicates
	echo " sambamba remove duplicates "
	singularity exec -B $mount_dir_singu $biotools sambamba markdup -t 8 -p --tmpdir $tmpdir --overflow-list-size 500000 $outdir/${idsample}.raw.sorted.bam $outdir/${idsample}.bam
	#7.
	for c in $(cat $chr_list); do
		echo " freebayes call variant per chromosome "
		singularity exec -B $mount_dir_singu $biotools freebayes -f $genome -r chr$c -g 2000  $outdir/${idsample}.bam | bgzip > $outdir/${idsample}.chr$c.vcf.gz
		tabix -p vcf $outdir/${idsample}.chr$c.vcf.gz
		zcat $outdir/${idsample}.chr$c.vcf.gz | singularity exec -B $mount_dir_singu $biotools vcffilter -f "QUAL > 20" > $outdir/${idsample}.chr$c.fb.filt.vcf
		bgzip $outdir/${idsample}.chr$c.fb.filt.vcf
		tabix -p vcf $outdir/${idsample}.chr$c.fb.filt.vcf.gz
		#10.
		echo " normalization "
		singularity exec -B $mount_dir_singu $biotools vt normalize -n $outdir/${idsample}.chr$c.fb.filt.vcf.gz -r $genome -o $outdir/${idsample}.chr$c.fb.norm.vcf.gz 
		tabix -p vcf $outdir/${idsample}.chr$c.fb.norm.vcf.gz
		#11.
		echo " decompose "
		singularity exec -B $mount_dir_singu $biotools vt decompose_blocksub  $outdir/${idsample}.chr$c.fb.norm.vcf.gz  -o $outdir/${idsample}.chr$c.fb.norm.decompose.vcf.gz
		tabix -p vcf $outdir/${idsample}.chr$c.fb.norm.decompose.vcf.gz

	;done

