#!/bin/bash
#from bam to fastq

# Positional arguments
sample_name=${@:$OPTIND:1}
file_name=${@:$OPTIND+1:1}
out_dir=${@:$OPTIND+2:1}

if test -f "$file_name"; then
	#sort
	echo "sorting file" $file_name
	samtools sort -n $file_name -o $out_dir/$sample_name.qsort.bam
	#convert
	echo "convert from sorted bam to fastq files R1 & R2" $file_name
	bedtools bamtofastq -i $out_dir/$sample_name.qsort.bam -fq $out_dir/$sample_name.R1.fq -fq2 $out_dir/$sample_name.R2.fq
	echo "deleting tmp file..." 
	rm $out_dir/$sample_name.qsort.bam
	echo "Done...exit" 
else 
	echo "ERROR: no input given...exit!"
	exit
fi