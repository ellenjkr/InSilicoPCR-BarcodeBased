#!/bin/bash

mkdir -p output

# for seq_file in input/*.gz; do  # Para cada arquivo de sequências
# 	primers=( input/*.tsv )
# 	while read name	forward	reverse min limit; do
# 		bash virtual_pcr.sh -c 0.9 -n $name -f $forward -r $reverse -i $seq_file -m $min -l $limit -t 4
# 	done < ${primers[0]};
# done



for primers in input/*.tsv; do
	seq_file=( input/*.gz )
	while IFS='	' read -r name	forward	reverse min limit; do
		bash virtual_pcr.sh -c "0.85" -n $name -f $forward -r $reverse -i $seq_file -m $min -M $limit -t 10;
	done < ${primers[0]};
done
