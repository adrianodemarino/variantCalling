#!/bin/bash
#https://bioinformaticsworkbook.org/dataAnalysis/VariantCalling/freebayes-dnaseq-workflow.html
#tutotial for variant calling
genome=/data/resources/genomes/human/hg38.fa
tmpdir=/scratch
biotools=/data/resources/img/biotools.img
cgt=/data/resources/img/CGT.img
mount_dir_singu=/data
outdir=/data/research/NGS/results/memorial_exome
sample_list=/data/research/NGS/results/memorial_exome/missing_sample.list



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
	freebayes -f $genome -r chr$c -g 2000  $outdir/$idsample/$idsample.bam | bgzip > $outdir/$idsample/$idsample.chr$c.vcf.gz && touch $idsample.chr$c_ok
	#tabix -p vcf $outdir/$idsample/$idsample.chr$c.vcf.gz
fi

echo "quality filter"
FILE=$outdir/$idsample/$idsample.chr$c.fb.filt.vcf
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	zcat $outdir/$idsample/$idsample.chr$c.vcf.gz | vcffilter -f "QUAL > 20" > $outdir/$idsample/$idsample.chr$c.fb.filt.vcf
	bgzip $outdir/$idsample/$idsample.chr$c.fb.filt.vcf
	#tabix -p vcf $outdir/$idsample/$idsample.chr$c.fb.filt.vcf.gz
fi

echo "normalize"
FILE=$outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	vt normalize -n $outdir/$idsample/$idsample.chr$c.fb.filt.vcf.gz -r $genome -o $outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz 
	#tabix -p vcf $outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz
fi

echo "decompose"
FILE=$outdir/$idsample/$idsample.chr$c.fb.norm.decompose.vcf.gz
if test -f "$FILE"; then
	echo ">>> $FILE exist"
else 
	vt decompose_blocksub  $outdir/$idsample/$idsample.chr$c.fb.norm.vcf.gz  -o $outdir/$idsample/$idsample.chr$c.fb.norm.decompose.vcf.gz
	#tabix -p vcf $outdir/$idsample/$idsample.chr$c.fb.norm.decompose.vcf.gz
fi




