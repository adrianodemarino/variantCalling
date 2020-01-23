#!/bin/bash
#launch Kore free-bayes in parallel
sample_list=/opt/NGS/data/priority_sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt

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
 192477)

for id in "${samples[@]}"; do
	echo "sample ID: $id"
	for c in $(cat $chr_list); do 
		/home/$(whoami)/github/variantCalling/kore-freebayes.sh -i $id -c $c & 
	done
done
