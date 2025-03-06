import sys
import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict


def adicionar_contador_ids(fasta_file):
    # Dicionário para contar os IDs repetidos
    id_count = defaultdict(int)

    # Ler o conteúdo do arquivo FASTA
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Lista para armazenar as sequências modificadas
    new_sequences = []

    for record in sequences:
        # Incrementar o contador de cada ID
        id_count[record.id] += 1

        # Se houver repetição, adicionar o contador ao ID
        if id_count[record.id] > 1:
            record.id = f"{record.id}_{id_count[record.id]}"
            record.description = record.id  # Atualizar a descrição com o novo ID

        # Adicionar o record modificado à lista de novas sequências
        new_sequences.append(record)

    # Sobrescrever o arquivo FASTA com os IDs modificados
    SeqIO.write(new_sequences, fasta_file, "fasta")

barcodes_file = sys.argv[1]
primer_name = sys.argv[2]
out_path = sys.argv[3]
primincl = sys.argv[4]
forwardp = sys.argv[5]
reversep = sys.argv[6]

# barcodes_file = 'barcodes/barcodes.tsv'
# primer_name = 'ANI2'
# out_path = 'output/lib8/test/primers_barcodes/'

df = pd.read_csv(barcodes_file, sep='\t', names=['SAMPLE_ID', 'PRIMER_ID', 'FORWARD_BARCODE', 'REVERSE_BARCODE'])
df = df.loc[df['PRIMER_ID'] == primer_name]

samples = {}
for sample in df['SAMPLE_ID']:
	samples[sample] = []

forward_sam = pd.read_csv('bbmap_out/forward.sam', sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])
forward_sam.loc[forward_sam['RNAME'] != '*', 'RNAME'] = forward_sam.loc[forward_sam['RNAME'] != '*', 'RNAME'].str.split().str[0]


reverse_sam = pd.read_csv('bbmap_out/reverse.sam', sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])
reverse_sam.loc[reverse_sam['RNAME'] != '*', 'RNAME'] = reverse_sam.loc[reverse_sam['RNAME'] != '*', 'RNAME'].str.split().str[0]

adicionar_contador_ids(f'{out_path}{primer_name}.fa')
record_dict = SeqIO.to_dict(SeqIO.parse(f'{out_path}{primer_name}.fa', "fasta"))
from Bio.SeqIO.QualityIO import FastqGeneralIterator
fastq_record = list(FastqGeneralIterator(open(f'{out_path}{primer_name}.fq')))
fastq_record_df = {'id': [], 'title': [], 'seq': [], 'qual': []}
for title, seq, qual in fastq_record:
    record_id = title.split(' ')[0]
    fastq_record_df['id'].append(record_id)
    fastq_record_df['title'].append(title)
    fastq_record_df['seq'].append(seq)
    fastq_record_df['qual'].append(qual)
fastq_record_df = pd.DataFrame(fastq_record_df)

data = {"ID": [seq_record.id for seq_record in record_dict.values()], "Description": [seq_record.description for seq_record in record_dict.values()], "Sequence": [str(seq_record.seq) for seq_record in record_dict.values()]}
sequences_df = pd.DataFrame(data)
sequences_df = sequences_df.merge(forward_sam, left_on='ID', right_on='RNAME')
sequences_df = sequences_df.merge(reverse_sam, left_on='ID', right_on='RNAME')

sequences_df = sequences_df[['ID', 'Sequence', 'QNAME_x', 'POS_x', 'QNAME_y', 'POS_y']]

if primincl == "f":
    sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'SEQ_F_BAR'] = sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[-5:].reverse_complement()))
    sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'SEQ_F_BAR'] = sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:5]))

    sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'] = sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:-5-len(forwardp)]))
    sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'] = sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[5+len(forwardp):]))

    sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'SEQ_R_BAR'] = sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[-5:].reverse_complement()))
    sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'SEQ_R_BAR'] = sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:5]))

    sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'] = sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:-5-len(reversep)].reverse_complement()))
    sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'] = sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[5+len(reversep):]))

else:
    sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'SEQ_F_BAR'] = sequences_df.loc[sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[-5:].reverse_complement()))
    sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'SEQ_F_BAR'] = sequences_df.loc[~sequences_df['QNAME_x'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:5]))

    sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'SEQ_R_BAR'] = sequences_df.loc[sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[-5:].reverse_complement()))
    sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'SEQ_R_BAR'] = sequences_df.loc[~sequences_df['QNAME_y'].str.contains('r_'), 'Sequence'].apply(lambda x: str(Seq(x)[:5]))


samples_sequences_df = sequences_df.merge(df, left_on=['SEQ_F_BAR', 'SEQ_R_BAR'], right_on=['FORWARD_BARCODE', 'REVERSE_BARCODE'])

grouped = samples_sequences_df.groupby('SAMPLE_ID')
for sample_id, data in grouped:
    file_name_fa = f"{out_path}{sample_id}.fa"
    file_name_fq = f"{out_path}{sample_id}.fastq.gz"


    with open(file_name_fa, 'w') as file:
        for index, row in data.iterrows():
            if 'r_' not in row['QNAME_x'] and 'r_' in row['QNAME_y'] and row['POS_x'] < row['POS_y']: 
                file.write(f">{row['SAMPLE_ID']}_seq{index}\n{row['Sequence']}\n")
            
            elif 'r_' not in row['QNAME_y'] and 'r_' in row['QNAME_x'] and row['POS_y'] < row['POS_x']: 
                # sequence = str(Seq(row['Sequence']).reverse_complement())
                sequence = row['Sequence']
                file.write(f">{row['SAMPLE_ID']}_seq{index}\n{sequence}\n")

    fastq_records = fastq_record_df.merge(data, left_on='id', right_on='ID').reset_index(drop=True)
    content = ''
    # with open(file_name_fq, 'w') as file:
    for index, row in fastq_records.iterrows():
        if 'r_' not in row['QNAME_x'] and 'r_' in row['QNAME_y'] and row['POS_x'] < row['POS_y']: 
            # file.write(f"@{row['title']}\n{row['seq']}\n+\n{row['qual']}\n")
            content += f"@{row['title']}_seq{index}\n{row['seq']}\n+\n{row['qual']}\n"
        
        elif 'r_' not in row['QNAME_y'] and 'r_' in row['QNAME_x'] and row['POS_y'] < row['POS_x']: 
            # sequence = str(Seq(row['seq']).reverse_complement())
            sequence = row['seq']
            # file.write(f"@{row['title']}\n{sequence}\n+\n{row['qual']}\n")
            content += f"@{row['title']}_seq{index}\n{sequence}\n+\n{row['qual']}\n"

    import gzip
    f = gzip.open(file_name_fq, 'wb')
    content = str.encode(content)
    f.write(content)
    f.close()







