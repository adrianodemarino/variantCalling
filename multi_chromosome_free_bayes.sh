#!/bin/bash
#launch Kore free-bayes in parallel
sample_list=/opt/NGS/data/sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt



for id in $(cat $sample_list); do 
	for c in $(cat $chr_list); do 
		/home/$(whoami)/github/variantCalling/kore-freebayes.sh -i $id -c $c & 
	done
done
