from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import  difflib


forward = sys.argv[1]
forward_comp = str(Seq(forward).reverse_complement())
reverse = sys.argv[2]
reverse_comp = str(Seq(reverse).reverse_complement())

sequences_fasta = sys.argv[3]

new_sequences = {'F_R': [], 'F_R_WRONG': [], 'FCR': [], 'RCR': []}
all_sequences = []

for pos, record in enumerate(SeqIO.parse(sequences_fasta, "fasta")):
	seq = record.seq
	
	left_primer = seq[:len(forward)]
	right_primer = seq[-len(reverse):]

	
	matches_left = difflib.get_close_matches(left_primer, [forward, forward_comp, reverse, reverse_comp])
	matches_right = difflib.get_close_matches(right_primer, [forward, forward_comp, reverse, reverse_comp])


	record = SeqRecord(
		record.seq,
		id=str(pos + 1),
		description=''
	)
	all_sequences.append(record)

	if matches_left[0] == forward:
		record = SeqRecord(
    		record.seq,
    		id=str(pos + 1),
    		description=''
		)
		new_sequences['F_R'].append(record)

	elif matches_left[0] == forward_comp:
		record = SeqRecord(
    		record.seq.complement(),
    		id=str(pos + 1),
    		description='C'
		)
		new_sequences['FCR'].append(record)

	elif matches_left[0] == reverse:
		record = SeqRecord(
    		record.seq.reverse_complement(),
    		id=str(pos + 1),
    		description='CR'
		)

		if matches_right[0] == forward:
			new_sequences['F_R_WRONG'].append(record)
		else:
			new_sequences['F_R'].append(record)

	elif matches_left[0] == reverse_comp:
		record = SeqRecord(
    		record.seq[::-1],
    		id=str(pos + 1),
    		description='R'
		)
		new_sequences['RCR'].append(record)
	
	else:
		print('Algo deu errado. Seq ID:', record.id)


values = [len(new_sequences[key]) for key in new_sequences.keys()]
total = sum(values)
percentages = [value * 100. / total for value in values]

print('Percentages:')
print(f'{round(percentages[0], 2)}% Forward or Reverse')
print(f'{round(percentages[1], 2)}% Forward - Reverse Complement')
print(f'{round(percentages[2], 2)}% Reverse - Reverse Complement')

for sequences_set in new_sequences.keys():
	file_name = sequences_fasta.split('/')[2]
	file_name = f'output/{sequences_set}/{sequences_set}_{file_name}'
	SeqIO.write(new_sequences[sequences_set], file_name, "fasta")


file_name = sequences_fasta.split('/')[2]
SeqIO.write(all_sequences, f'output/All/{file_name}', "fasta")
