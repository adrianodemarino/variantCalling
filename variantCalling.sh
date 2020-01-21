#!/bin/bash
#pipe for varinat calling
genome=/opt/NGS/genomes/hg38.fa.gz
genome2=/opt/NGS/genomes/hg38.p12.fa
tmpdir=/opt/NGS/scratch/
biotools=/home/adriano/src/img/biotools.img
mount_dir_singu=/opt/NGS/
inputdir=/opt/NGS/data/memorial_hospital
outdir=/opt/NGS/results/memorial_hospital
sample_list=/opt/NGS/data/sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt

for idsample in $(cat $sample_list); do 
	mkdir $outdir/$idsample
	cd $outdir/$idsample
	idsample_R1=`ls $inputdir/*$idsample*_1*`
	idsample_R2=`ls $inputdir/*$idsample*_2*`
	echo Start analysis of sample $idsample_R1 $idsample_R2
	#1. Dowload human hg38 genome
	#2. index it with: `bwa index -a bwtsw hg38.fa.gz`
	#3. mapping
	echo " mapping bwa "
	bwa mem -t 8 -R "@RG\tID:$idsample\tSM:$idsample" $genome $idsample_R1 $idsample_R2 | singularity exec -B $mount_dir_singu $biotools samtools view -b - > $outdir/$idsample.raw.bam
	#4. View
	#singularity exec -B $mount_dir_singu $biotools samtools view -h $outdir/$idsample.raw.bam | less -S
	#5. 
	echo " sambamba sort "
	singularity exec -B $mount_dir_singu $biotools sambamba sort -t 10 -m 25G --tmpdir $tmpdir -o $outdir/${idsample}.raw.sorted.bam $outdir/${idsample}.raw.bam
	#6. remove duplicates
	echo " sambamba remove duplicates "
	singularity exec -B $mount_dir_singu $biotools sambamba markdup -t 8 -p --tmpdir $tmpdir --overflow-list-size 500000 $outdir/${idsample}.raw.sorted.bam $outdir/${idsample}.bam
	#7.

done

