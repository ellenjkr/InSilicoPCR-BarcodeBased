import sys
import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
# barcodes_file = sys.argv[1]
# primer_name = sys.argv[2]

barcodes_file = 'barcodes/barcodes.tsv'
primer_name = 'BAC1'


df = pd.read_csv(barcodes_file, sep='\t', names=['SAMPLE_ID', 'PRIMER_ID', 'FORWARD_BARCODE', 'REVERSE_BARCODE'])
df = df.loc[df['PRIMER_ID'] == primer_name]

samples = {}
for sample in df['SAMPLE_ID']:
	samples[sample] = []

forward_sam = pd.read_csv('bbmap_out/forward.sam', sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])
forward_sam.loc[forward_sam['RNAME'] != '*', 'RNAME'] = forward_sam.loc[forward_sam['RNAME'] != '*', 'RNAME'].str.split().str[0]


reverse_sam = pd.read_csv('bbmap_out/reverse.sam', sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])
reverse_sam.loc[reverse_sam['RNAME'] != '*', 'RNAME'] = reverse_sam.loc[reverse_sam['RNAME'] != '*', 'RNAME'].str.split().str[0]

record_dict = SeqIO.to_dict(SeqIO.parse(f'output/{primer_name}.fa', "fasta"))


data = {"ID": [seq_record.id for seq_record in record_dict.values()], "Description": [seq_record.description for seq_record in record_dict.values()], "Sequence": [str(seq_record.seq) for seq_record in record_dict.values()]}
sequences_df = pd.DataFrame(data)
sequences_df = sequences_df.merge(forward_sam, left_on='ID', right_on='RNAME')
sequences_df = sequences_df.merge(reverse_sam, left_on='ID', right_on='RNAME')

sequences_df = sequences_df[['ID', 'Sequence', 'QNAME_x', 'POS_x', 'QNAME_y', 'POS_y']]

sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'SEQ_F_BAR'] = sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[-5:].reverse_complement()))
sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'SEQ_F_BAR'] = sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:5]))

sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'SEQ_R_BAR'] = sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[-5:].reverse_complement()))
sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'SEQ_R_BAR'] = sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:5]))

samples_sequences_df = sequences_df.merge(df, left_on=['SEQ_F_BAR', 'SEQ_R_BAR'], right_on=['FORWARD_BARCODE', 'REVERSE_BARCODE'])

grouped = samples_sequences_df.groupby('SAMPLE_ID')
for sample_id, data in grouped:
    file_name = f"output/{sample_id}.fa"
    with open(file_name, 'w') as file:
        for index, row in data.iterrows():

            if 'r_' not in row['QNAME_x'] and 'r_' in row['QNAME_y'] and row['POS_x'] < row['POS_y']: 
                file.write(f">{row['ID']}\n{row['Sequence']}\n")
            
            elif 'r_' not in row['QNAME_y'] and 'r_' in row['QNAME_x'] and row['POS_y'] < row['POS_x']: 
                sequence = str(Seq(row['Sequence']).reverse_complement())
                file.write(f">{row['ID']}\n{sequence}\n")

