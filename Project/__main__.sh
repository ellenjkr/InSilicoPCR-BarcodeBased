#!/bin/bash

mkdir -p output
mkdir -p output/All
mkdir -p output/F_R
mkdir -p output/F_R_WRONG
mkdir -p output/FCR
mkdir -p output/RCR

for seq_file in input/*.gz; do  # Para cada arquivo de sequências
	primers=( input/*.tsv )
	while read name	forward	reverse min limit; do
		bash virtual_pcr.sh -c 0.9 -n $name -f $forward -r $reverse -i $seq_file -m $min -l $limit
	done < ${primers[0]};
done
