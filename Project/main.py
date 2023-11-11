import yaml
import os
import pandas as pd
import subprocess
import boxplot


def read_config_file(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


config = read_config_file('config.yaml')

if config["HAS_BARCODES"] is True:
    config["HAS_BARCODES"] = 'true'
elif config["HAS_BARCODES"] is False:
    config["HAS_BARCODES"] = 'false'

if config['OUTPUT_PATH'][-1] != '/':
    config['OUTPUT_PATH'] += '/'


if not os.path.exists(config['OUTPUT_PATH']): 
    os.makedirs(config['OUTPUT_PATH']) 

input_path = config['INPUT_PATH'] + '/' if config['INPUT_PATH'][-1] != '/' else config['INPUT_PATH']

primers_files = []
seq_file = None
for file in os.listdir(config['INPUT_PATH']):
    if '.tsv' in file:
        primers_files.append(input_path + file)
    elif '.gz' in file or '.fq' in file or '.fastq' in file:
        seq_file = input_path + file

primer_names = []
for primer_file in primers_files:
    primer_set = primer_file.split('.')[0]
    primer_set = primer_set.split('/')[-1]
    
    df = pd.read_csv(primer_file, sep='\t', names=['name', 'forward', 'reverse', 'min', 'max'])
    for idx, row in df.iterrows():
        if row['name'] not in primer_names:
            primer_names.append(row['name'])
        subprocess.run(
            f'bash virtual_pcr.sh -s {primer_set} -o {config["OUTPUT_PATH"]} -c "{config["IDENTITY"]}" -b {config["HAS_BARCODES"]} -l {config["BARCODES_LEN"]} -B {config["BARCODES_PATH"]} -n {row["name"]} -f {row["forward"]} -r {row["reverse"]} -i {seq_file} -m {row["min"]} -M {row["max"]} -t {config["THREADS"]}',
            shell=True,
            executable='/bin/bash'
        )

boxplot.run(config, primer_names)
