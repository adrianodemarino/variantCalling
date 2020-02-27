#!/bin/bash
#launch Kore free-bayes in parallel
sample_list=/data/research/NGS/results/memorial_exome/missing_sample.list


for id in `cat $sample_list`; do
    echo "sample ID: $id"
    for c in {1..22} X Y; do 
        ~/repositories/illumina-pipeline/kore-freebayes.sh -i $id -c $c &
    done
    wait
done
