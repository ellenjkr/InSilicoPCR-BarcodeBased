import sys
import pandas as pd


def change_position(sam_file, barcode_length):
	alignments_df = pd.read_csv(sam_file, sep='\t', skiprows=[0], names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'Identity'])
	mask = (~alignments_df['QNAME'].str.contains('r_')) & (alignments_df['RNAME'] != '*')
	alignments_df.loc[mask, 'POS'] = alignments_df.loc[mask, 'POS'] - barcode_length

	mask = (alignments_df['QNAME'].str.contains('r_')) & (alignments_df['RNAME'] != '*')
	alignments_df.loc[mask, 'POS'] = alignments_df.loc[mask, 'POS'] + barcode_length

	alignments_df.to_csv(sam_file, sep='\t', header=False, index=False)

	with open(sam_file, 'r') as sam_temp:
		sam_content = sam_temp.read()

	with open(sam_file, 'w') as sam_temp:
		sam_temp.write('@HD	VN:1.4	SO:unsorted\n')
		sam_temp.write(sam_content)


change_position(sys.argv[1], 5)
change_position(sys.argv[2], 5)


