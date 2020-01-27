#!/bin/bash
#launch Kore free-bayes in parallel
sample_list=/opt/NGS/data/priority_sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt

samples=(
192361
192824)


for id in "${samples[@]}"; do
    echo "sample ID: $id"
    for c in $(cat $chr_list); do 
        /home/$(whoami)/github/variantCalling/kore-freebayes.sh -i $id -c $c &
    done
    wait
done
