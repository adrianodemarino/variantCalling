#!/bin/bash
#from bam to fastq

# Positional arguments
sample_name=${@:$OPTIND:1}
file_name=${@:$OPTIND+1:1}
out_dir=${@:$OPTIND+2:1}

#sort file
samtools sort -n $file_name -o $out_dir/$sample_name.qsort.bam
#convert
bedtools bamtofastq -i $out_dir/$sample_name.qsort.bam -fq $out_dir/$sample_name.r1.fq -fq2 $out_dir/$sample_name.r2.fq
#remove tmp file
rm $out_dir/$sample_name.qsort.bam


