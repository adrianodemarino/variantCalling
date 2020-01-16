#pipe for varinat calling
idsample=m192244-19WES
1. Dowload human hg38 genome
2. index it with: `bwa index -a bwtsw hg38.fa.gz`
3. mapping: bwa mem -t 8 -R "@RG\tID:$idsample\tSM:$idsample" /home/adriano/Documents/data/genomes/hg38.p12.fa m192244-19WES_L4_1.fq.gz m192244-19WES_L4_2.fq.gz | samtools view -b - > $idsample.raw.bam
4. singularity exec -B /home/adriano ~/src/img/biotools.img samtools view -h m192244-19WES.raw.bam | less -S
5. singularity exec -B /home/adriano ~/src/img/biotools.img sambamba sort -t 10 -m 20G --tmpdir /home/adriano/Documents/data/seq_memorial/scratch -o ${idsample}.raw.sorted.bam ${idsample}.raw.bam
6. singularity exec -B /home/adriano ~/src/img/biotools.img sambamba markdup -t 8 -p --tmpdir scratch --overflow-list-size 500000 ${idsample}.raw.sorted.bam ${idsample}.bam
7. singularity exec -B /home/adriano ~/src/img/biotools.img freebayes -f /home/adriano/Documents/data/genomes/hg38.p12.fa -g 2000  ${idsample}.bam | bgzip > ${idsample}.vcf.gz
8. tabix -p vcf ${idsample}.vcf.gz 
9. zcat ${idsample}.vcf.gz | singularity exec -B /home/adriano ~/src/img/biotools.img vcffilter -f "QUAL > 20" > ${idsample}.fb.filt.vcf | bgzip ${idsample}.fb.filt.vcf
10. singularity exec -B /home/adriano ~/src/img/biotools.img vt normalize ${idsample}.fb.filt.vcf.gz -r /home/adriano/Documents/data/genomes/hg38.p12.fa -o  ${idsample}.fb.norm.vcf.gz
11. singularity exec -B /home/adriano ~/src/img/biotools.img vt decompose_blocksub  ${idsample}.fb.norm.vcf.gz  -o ${idsample}.fb.norm.decompose.vcf.gz
