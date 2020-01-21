#pipe for varinat calling
#!/bin/bash
genome=/opt/NGS/genomes/hg38.p12.fa
tmpdir=/opt/NGS/scratch
biotools=/home/adriano/src/img/biotools.img
mount_dir_singu=/home/adriano
for word in $(cat sample_list.txt); 
	
	do 
	idsample_R1=`ls *$word*_1*`
	idsample_R2=`ls *$word*_2*`
	echo Start analysis of sample $idsample_R1 $idsample_R2
	#1. Dowload human hg38 genome
	#2. index it with: `bwa index -a bwtsw hg38.fa.gz`
	#3. mapping
	bwa mem -t 8 -R "@RG\tID:$idsample\tSM:$idsample"  $idsample_R1 $idsample_R2 | samtools view -b - > $idsample.raw.bam
	#4.
	singularity exec -B $mount_dir_singu $biotools samtools view -h m192244-19WES.raw.bam | less -S
	#5.
	singularity exec -B $mount_dir_singu $biotools sambamba sort -t 10 -m 20G --tmpdir $tmpdir -o ${idsample}.raw.sorted.bam ${idsample}.raw.bam
	#6.
	singularity exec -B $mount_dir_singu $biotools sambamba markdup -t 8 -p --tmpdir scratch --overflow-list-size 500000 ${idsample}.raw.sorted.bam ${idsample}.bam
	#7.
	singularity exec -B $mount_dir_singu $biotools freebayes -f $genome -g 2000  ${idsample}.bam | bgzip > ${idsample}.vcf.gz
	#8.
	tabix -p vcf ${idsample}.vcf.gz 
	#9.
	zcat ${idsample}.vcf.gz | singularity exec -B $mount_dir_singu $biotools vcffilter -f "QUAL > 20" > ${idsample}.fb.filt.vcf
	bgzip ${idsample}.fb.filt.vcf
	tabix -p vcf ${idsample}.fb.filt.vcf.gz
	#10.
	singularity exec -B $mount_dir_singu $biotools vt normalize -n ${idsample}.fb.filt.vcf.gz -r $genome -o ${idsample}.fb.norm.vcf.gz 
	tabix -p vcf ${idsample}.fb.norm.vcf.gz
	#11.
	singularity exec -B $mount_dir_singu $biotools vt decompose_blocksub  ${idsample}.fb.norm.vcf.gz  -o ${idsample}.fb.norm.decompose.vcf.gz
	tabix -p vcf ${idsample}.fb.norm.decompose.vcf.gz
	
	;done

