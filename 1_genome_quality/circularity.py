import os, Bio.SeqIO, sys, argparse
from Bio.Seq import Seq
import string

tab = string.maketrans("ACTG", "TGAC")

def reverse_complement(seq):
    return seq.translate(tab)[::-1]

def compute_fwd_overlap(seq, revcomp=False):
	max_overlap = min(10000, len(seq)/2)
	target = str(Seq(seq).reverse_complement()) if revcomp else seq
	for i in range(1,max_overlap)[::-1]:
		if seq[0:i] == target[-i:]:
			return i
	return 0

def compute_rev_overlap(seq):
	max_overlap = min(10000, len(seq)/2)
	for i in range(1,max_overlap)[::-1]:
		if seq[0:i] == reverse_complement(seq[-i:]):
			return i
	return 0

def main(in_fna):
	data = {}
	handle = Bio.SeqIO.parse(in_fna, 'fasta')
	for index, seq in enumerate(handle):
		repeat = compute_fwd_overlap(str(seq.seq).upper())
		data[seq.id] = repeat
	return data


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('--in_fna', required=True, metavar='PATH')
	parser.add_argument('--out_tsv', required=True, metavar='PATH')
	parser.add_argument('--repeat_length', required=True, metavar='INT', type=int)
	args = vars(parser.parse_args())
	
	data = main(args['in_fna'])
	out = open(args['out_tsv'], 'w')
	out.write('contig_id\trepeat_length\tis_circular\n')
	for id, repeat in data.items():
		is_circular = 'Yes' if repeat >= args['repeat_length'] else 'No'
		row = [id, str(repeat), is_circular]
		out.write('\t'.join(row)+'\n')




