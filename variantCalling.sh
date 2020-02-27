#!/bin/bash
#pipe for varinat calling
genome=/data/resources/genomes/human/hg38.fa
tmpdir=/scratch
biotools=/home/adriano/img/biotools.img
mount_dir_singu=/data
inputdir=/data/research/NGS/data/clinical_exome_rappresentative/fastq
outdir=/data/research/NGS/data/clinical_exome_rappresentative/coverage
sample_list=/data/research/NGS/results/memorial_exome/missing_sample.list

for idsample in `cat $sample_list`; do
    dir=$outdir/$idsample
    if [ -d "$dir" ]; then
        echo ">>> $dir exist"
        echo ">>> SKIP"
    else 
        mkdir $outdir/$idsample
        idsample_R1=`ls $inputdir/*$idsample* | head -1`
        idsample_R2=`ls $inputdir/*$idsample* | tail -1`
        echo Start analysis of sample $idsample_R1 $idsample_R2
        #1. Dowload human hg38 genome
        #2. index it with: `bwa index -a bwtsw hg38.fa.gz`
        #3. mapping
        echo "> mapping bwa "
        bwa mem -t 30 -R "@RG\tID:$idsample\tSM:$idsample" $genome $idsample_R1 $idsample_R2 |  samtools view -b - > $outdir/$idsample/$idsample.raw.bam
        #4. View
        #singularity exec -B $mount_dir_singu $biotools samtools view -h $outdir/$idsample/$idsample.raw.bam | less -S
        #5. 
        echo "> sambamba sort "
        sambamba sort -t 30 -m 25G --tmpdir $tmpdir -o $outdir/$idsample/${idsample}.raw.sorted.bam $outdir/$idsample/${idsample}.raw.bam
        #6. mark duplicates
        echo "> sambamba remove duplicates "
        sambamba markdup -t 30 -p --tmpdir $tmpdir --overflow-list-size 500000 $outdir/$idsample/${idsample}.raw.sorted.bam $outdir/$idsample/${idsample}.bam
        #7.
    fi
done