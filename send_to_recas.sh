#!/bin/bash
#send file to recas

cd /media/adriano/TOSHIBA/memorial_exome_results/results/memorial_hospital

samples=( 182088
 190116
 190742
 190899
 191397
 192130
 192186
 192486
 192805
 192897
 182151
 190001
 190748
 190923
 191460
 192152
 192502
 192991
 182189
 190169
 190798
 191093
 191725
 192167
 192311
 193012
 182204
 190210
 190810
 191386
 192128
 192184
 192477
 192695)

for idsample in "${samples[@]}"; do 
	scp -r $idsample/$idsample*chr*deco*.vcf.gz enza@frontend.recas.ba.infn.it:/lustre/home/enza/adriano/project/memorial_hospital/results/variantcalling ; done
