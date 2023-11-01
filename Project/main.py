import yaml
import os
import pandas as pd
import subprocess


def read_config_file(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


config = read_config_file('config.yaml')

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

for primer_file in primers_files:
    df = pd.read_csv(primer_file, sep='\t', names=['name', 'forward', 'reverse', 'min', 'max'])
    for idx, row in df.iterrows():
        subprocess.run(
            f'bash virtual_pcr.sh -c "{config["IDENTITY"]}" -b {config["BARCODES_LEN"]} -n {row["name"]} -f {row["forward"]} -r {row["reverse"]} -i {seq_file} -m {row["min"]} -M {row["max"]} -t {config["THREADS"]}',
            shell=True,
            executable='/bin/bash'
        )
