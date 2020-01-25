#!/bin/bash
#send file to recas

cd /media/adriano/TOSHIBA/memorial_exome_results/results/memorial_hospital
#male infertility samples
samples=( 182088
182189
182204
190001
190169
190748
190899
191093
191725
192167
192184
192311
192477
192486
192502
192991)

for idsample in "${samples[@]}"; do 
	scp -r $idsample/$idsample*chr*deco*.vcf.gz enza@frontend.recas.ba.infn.it:/lustre/home/enza/adriano/project/memorial_hospital/results/variantcalling ; done

