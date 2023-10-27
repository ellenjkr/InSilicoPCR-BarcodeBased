from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import  difflib
import argparse


parser=argparse.ArgumentParser()

parser.add_argument('-f', type=str)
parser.add_argument('-r', type=str)
parser.add_argument('-s', type=str)
parser.add_argument('--f_r', action='store_true')
parser.add_argument('--r_fc', action='store_true')
parser.add_argument('--f_rc', action='store_true')
parser.add_argument('--fc_r', action='store_true')

args = parser.parse_args()


forward = args.f
forward_comp = str(Seq(forward).reverse_complement())
reverse = args.r
reverse_comp = str(Seq(reverse).reverse_complement())

sequences_fasta = args.s

f_r_orientation = args.f_r
r_fc_orientation = args.r_fc
f_rc_orientation = args.f_rc
fc_r_orientation = args.fc_r


sequences2save = []
for pos, record in enumerate(SeqIO.parse(sequences_fasta, "fasta")):
	seq = record.seq
	
	# forward -> reverse
	left_primer_fr = seq[:len(forward)]
	right_primer_fr = seq[-len(reverse):]

	
	matches_left_fr = difflib.get_close_matches(left_primer_fr, [forward, forward_comp, reverse, reverse_comp])
	best_match_left_fr = matches_left_fr[0]

	matches_right_fr = difflib.get_close_matches(right_primer_fr, [forward, forward_comp, reverse, reverse_comp])
	best_match_right_fr = matches_right_fr[0]

	# reverse -> forward
	left_primer_rf = seq[:len(reverse)]
	right_primer_rf  = seq[-len(forward):]

	
	matches_left_rf = difflib.get_close_matches(left_primer_rf, [forward, forward_comp, reverse, reverse_comp])
	best_match_left_rf = matches_left_rf[0]

	matches_right_rf = difflib.get_close_matches(right_primer_rf, [forward, forward_comp, reverse, reverse_comp])
	best_match_right_rf = matches_right_rf[0]


	# comparison of different orientations matches
	left_fr_score = difflib.SequenceMatcher(None, left_primer_fr, best_match_left_fr).ratio()
	right_fr_score = difflib.SequenceMatcher(None, right_primer_fr, best_match_right_fr).ratio()
	left_rf_score = difflib.SequenceMatcher(None, left_primer_rf, best_match_left_rf).ratio()
	right_rf_score = difflib.SequenceMatcher(None, right_primer_rf, best_match_right_rf).ratio()


	left_match = best_match_left_fr if left_fr_score > left_rf_score else best_match_left_rf
	right_match = best_match_right_fr if right_fr_score > right_rf_score else best_match_right_rf
	
	# check orientations

	record = SeqRecord(
		record.seq,
		id=str(pos + 1),
		description=''
	)


	if left_match == forward and right_match == reverse and f_r_orientation is True:
		record = SeqRecord(
    		record.seq,
    		id=str(pos + 1),
    		description=''
		)
		sequences2save.append(record)

	elif left_match == reverse and right_match == forward_comp and r_fc_orientation is True:
		record = SeqRecord(
    		record.seq,
    		id=str(pos + 1),
    		description=''
		)
		sequences2save.append(record)

	elif left_match == forward and right_match == reverse_comp and f_rc_orientation is True:
		record = SeqRecord(
    		record.seq,
    		id=str(pos + 1),
    		description=''
		)
		sequences2save.append(record)

	elif left_match == forward_comp and right_match == reverse and fc_r_orientation is True:
		record = SeqRecord(
    		record.seq,
    		id=str(pos + 1),
    		description=''
		)
		sequences2save.append(record)

SeqIO.write(sequences2save, sequences_fasta, "fasta")

