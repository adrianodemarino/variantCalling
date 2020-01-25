#!/bin/bash
#launch Kore free-bayes in parallel
sample_list=/opt/NGS/data/priority_sample_list.txt
chr_list=/opt/NGS/data/chr_list.txt

samples=(
181755
182151
190210
190423
190810
191080
191386
191397
192127
192130
192152
192186
192361
192369
192682
192805
192824
192897
193012)


for id in "${samples[@]}"; do
    echo "sample ID: $id"
    for c in $(cat $chr_list); do 
        /home/$(whoami)/github/variantCalling/kore-freebayes.sh -i $id -c $c &
    done
    wait
done
