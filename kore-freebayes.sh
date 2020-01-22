#!/bin/bash
#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/freebayes-dnaseq-workflow.html
#tutotial for variant calling
genome=/opt/NGS/genomes/hg38.fa.gz
genome2=/opt/NGS/genomes/hg38.p12.fa
tmpdir=/opt/NGS/scratch/
biotools=/home/adriano/src/img/biotools.img
mount_dir_singu=/opt/NGS/
inputdir=/opt/NGS/data/memorial_hospital
outdir=/opt/NGS/results/memorial_hospital
sample_list=/opt/NGS/data/priority_sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt

while getopts i:c: option
do
case "${option}"
in
i) idsample=${OPTARG};; #idsample
c) c=${OPTARG};; #chromosome
esac
done


echo "freebayes"

FILE=$outdir/$idsample/$idsample.chr$c.vcf.gz
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	singularity exec -B $mount_dir_singu $biotools freebayes -f $genome2 -r chr$c -g 2000  $outdir/$idsample/$idsample.bam | bgzip > $outdir/$idsample/$idsample.chr$c.vcf.gz
	tabix -p vcf $outdir/$idsample/$idsample.chr$c.vcf.gz
fi

echo "quality filter"
FILE=$outdir/$idsample/$idsample.chr$c.fb.filt.vcf
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	zcat $outdir/$idsample/$idsample.chr$c.vcf.gz | singularity exec -B $mount_dir_singu $biotools vcffilter -f "QUAL > 20" > $outdir/$idsample/$idsample.chr$c.fb.filt.vcf
	bgzip $outdir/$idsample/$idsample.chr$c.fb.filt.vcf
	tabix -p vcf $outdir/$idsample/$idsample.chr$c.fb.filt.vcf.gz
fi

echo "normalize"
FILE=$outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	singularity exec -B $mount_dir_singu $biotools vt normalize -n $outdir/$idsample/$idsample.chr$c.fb.filt.vcf.gz -r $genome2 -o $outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz 
	tabix -p vcf $outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz
fi

echo "decompose"
FILE=$outdir/$idsample/$idsample.chr$c.fb.norm.decompose.vcf.gz
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	singularity exec -B $mount_dir_singu $biotools vt decompose_blocksub  $outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz  -o $outdir/$idsample/$idsample.chr$c.fb.norm.decompose.vcf.gz
	tabix -p vcf $outdir/$idsample/$idsample.chr$c.fb.norm.decompose.vcf.gz
fi




