#!/bin/bash
#pipe for varinat calling
genome=/media/adriano/TOSHIBA/memorial_exome_results/genomes/hg38.fa.gz
genome2=/media/adriano/TOSHIBA/memorial_exome_results/genomes/hg38.p12.fa
tmpdir=/media/adriano/TOSHIBA/memorial_exome_results/scratch
biotools=/home/adriano/src/img/biotools.img
mount_dir_singu=/media/adriano/TOSHIBA/
inputdir=/opt/NGS/data/memorial_hospital
outdir=/media/adriano/TOSHIBA/memorial_exome_results/results/memorial_hospital
sample_list=/opt/NGS/data/priority_sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt

for idsample in 192361 192824; do
    dir=$outdir/$idsample
    if [ -d "$dir" ]; then
        echo ">>> $dir exist"
        echo ">>> SKIP"
    else 
        mkdir $outdir/$idsample
        idsample_R1=`ls $inputdir/*$idsample*_1*`
        idsample_R2=`ls $inputdir/*$idsample*_2*`
        echo Start analysis of sample $idsample_R1 $idsample_R2
        #1. Dowload human hg38 genome
        #2. index it with: `bwa index -a bwtsw hg38.fa.gz`
        #3. mapping
        echo "> mapping bwa "
        bwa mem -t 8 -R "@RG\tID:$idsample\tSM:$idsample" $genome $idsample_R1 $idsample_R2 | singularity exec -B $mount_dir_singu $biotools samtools view -b - > $outdir/$idsample/$idsample.raw.bam
        #4. View
        #singularity exec -B $mount_dir_singu $biotools samtools view -h $outdir/$idsample/$idsample.raw.bam | less -S
        #5. 
        echo "> sambamba sort "
        singularity exec -B $mount_dir_singu $biotools sambamba sort -t 10 -m 25G --tmpdir $tmpdir -o $outdir/$idsample/${idsample}.raw.sorted.bam $outdir/$idsample/${idsample}.raw.bam
        #6. mark duplicates
        echo "> sambamba remove duplicates "
        singularity exec -B $mount_dir_singu $biotools sambamba markdup -t 8 -p --tmpdir $tmpdir --overflow-list-size 500000 $outdir/$idsample/${idsample}.raw.sorted.bam $outdir/$idsample/${idsample}.bam
        #7.
    fi
done