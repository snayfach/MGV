#!/usr/bin/env python

import sys, os, time, shutil, sys, re
from collections import Counter
import numpy as np
import subprocess as sp

def run_pilercr(fasta, tmpdir):
	cmd = "bin/pilercr -in %s -out /dev/stdout" % fasta
	p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	o, e = p.communicate()
	c = p.returncode
	return str(o,encoding='utf-8'), e, c

def string_to_array(string):
	
	# initialize new array
	array = Array()
	array.contig_id = string.split('\n')[1].lstrip('>')
	array.reference_repeat = string.rstrip('\n').split('\n')[-1].split()[-1].upper()
	
	# find where repeat and spacers start in string
	header = re.split("\n=.*\n", string)[0].rstrip('\n').split('\n')[-1]
	repeat_string_index = header.rfind("Repeat")
	spacer_string_index = header.rfind("Spacer")
	
	# split up string per repeat/spacer pair
	lines = re.split("\n=.*\n", string)[1].rstrip('\n').split('\n')
	for line_index, line in enumerate(lines):
	
		# initialize new array element
		array_item = ArrayElement()
		
		# get individual spacers
		if len(line) < spacer_string_index + 1: # string end before spacer starts
			array_item.spacer_seq = ''
		else: # spacer sequence present
			array_item.spacer_seq = line[spacer_string_index:].split()[0].upper()
		
		# get individual repeats
		array_item.repeat_start = int(line.split()[0])
		array_item.repeat_seq = ''
		array_item.gaps = 0
		repeat_string = line[repeat_string_index:].split()[0].upper()
		for i, c in enumerate(repeat_string):
			if c == '.': # reference
				array_item.repeat_seq += array.reference_repeat[i]
			elif c == '-': # indel
				array_item.gaps += 1
				array_item.repeat_seq += c
				continue
			else: # substitution
				array_item.repeat_seq += c
		
		# repeat length does not include gaps
		array_item.repeat_length = len(array_item.repeat_seq.replace('-', ''))
		
		# store array element
		array.items.append(array_item)
	
	# last spacer is always empty
	array.items[-1].spacer_seq = ''
	return array

def parse_pilercr(out, contigs):
	
	arrays = []
	if out.find('DETAIL REPORT\n') == -1:
		return arrays

	strings = out.split('DETAIL REPORT\n')[1].split('SUMMARY BY SIMILARITY\n')[0].split('\n\n\nArray ')[1:]

	for string in strings:
		
		# parse array from output
		array = string_to_array(string)
		array.id = len(arrays)
		
		# define consnsus repeat
		array.define_consensus()
		
		# fix repeat start positions; pilercr counts gaps towards repeat starts
		sum_gaps = 0
		for i in range(1, len(array.items)):
			sum_gaps += array.items[i-1].gaps
			array.items[i].repeat_start -= sum_gaps
		
		arrays.append(array)

	return arrays

def run_crt(fasta, tmpdir, xmx):
	""" Run CRT; specify mem usage; return output; exit if non-zero exit code
	"""
	cmd = "java -Xmx%sm -cp bin/CRT1.2-CLI.jar crt %s /dev/stdout" % (xmx, fasta)
	p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	o, e = p.communicate()
	c = p.returncode
	return str(o,encoding='utf-8'), e, c

def parse_crt(out, contigs):
	""" Parse CRT output; make sure arrays do not span 2 contigs
	"""
	arrays = []
	contig_index = 0
	
	# No CRISPRs; return empty list
	if out.find('\nCRISPR') == -1:
		return arrays

	# Each string is a CRISPR array
	for string in out.split('\nCRISPR ')[1:]:
		
		array = Array()
		array.id = len(arrays)
		
		for line in string.split('\n')[3:]:
			
			row = line.rstrip().split()
			if line[0] == '-': break # end of array
			
			array_item = ArrayElement()
			
			# parse some element attributes from output
			array_item.repeat_start = int(row[0])
			array_item.repeat_seq = row[1].upper()
			array_item.repeat_length = len(array_item.repeat_seq.replace('-', ''))
			array_item.spacer_seq = row[2].upper() if len(row) > 2 else ""
			
			# determine element coordinates:
			# a) contig based on the repeat start position
			while True:
				if array_item.repeat_start > contigs[contig_index].end:
					contig_index += 1
				else:
					array_item.contig_id = contigs[contig_index].id
					break
			# b) repeat/spacer start positions relative to contig
			array_item.repeat_start = array_item.repeat_start - contigs[contig_index].start + 1
			array_item.spacer_start = array_item.repeat_start + len(array_item.repeat_seq)

			# set contig id for Array
			# only do this if repeat doesnt start or end with N
			# if it does start with N; then the start of repeat may come from upstream contig
			if array_item.repeat_seq[0] != 'N' and array_item.repeat_seq[-1] != 'N':
				array.contig_id = contigs[contig_index].id
				array.contig_index = contig_index
			array.items.append(array_item)

		# in CRT output, some arrays span >1 contig;
		# in these cases, delete repeat/spacer pairs that overlap 1st contig
		del_items = []
		for item_index, item in enumerate(array.items):
			# repeat starts on previous contig
			if (item.contig_id != array.contig_id or
					# repeat end on next contig
					item.repeat_start + len(item.repeat_seq) >= contigs[array.contig_index].length or
					# spacer ends on next contig
					item.spacer_start + len(item.spacer_seq) >= contigs[array.contig_index].length
				):
				del_items.append(item_index)
		for item_index in del_items[::-1]:
			del array.items[item_index]

		# define consensus repeat
		array.define_consensus()

		arrays.append(array)

	return arrays

class Array:
	def __init__(self):
		self.contig = None
		self.start = None
		self.end = None
		self.items = []
		self.reference_repeat = None
		self.consensus_repeat = None

	def prune(self, max_nfraction=0.80, warn=False):
		""" Delete array elements where repeat contains too many Ns
			to do: only apply to leading/trailing?
		"""
		indexes = []
		for index, item in enumerate(self.items):
			repeat_n = item.repeat_seq.count('N')
			nfraction = repeat_n/float(item.repeat_length)
			if nfraction > max_nfraction:
				indexes.append(index)
		if len(indexes) > 0 and warn:
			sys.stderr.write("Warning: Too many Ns for array\n")
		for index in indexes[::-1]:
			del self.items[index]
					
	def replace_leading_ns(self, append_ns, contig, warn=True, replace_char='_'):
		""" Replaces "leftmost overhanging Ns" in leading repeat with - character
		    No need to replace Ns in spacer: would mean that entire repeat was Ns and these elements are dropped
		"""
		self.items[0].left_overhang = max(append_ns - self.items[0].repeat_start + 1, 0)
		if self.items[0].left_overhang > 0:
			if warn: 
				sys.stderr.write("Warning: Leading Ns for array\n")
				sys.stderr.write("  before:%s\n" % self.items[0].repeat_seq)			
			seq = list(self.items[0].repeat_seq)			
			for i in range(self.items[0].left_overhang):
				seq[i] = replace_char
			self.items[0].repeat_seq = ''.join(seq)
			if warn: 
				sys.stderr.write("   after:%s\n\n" % self.items[0].repeat_seq)			
			
	def replace_trailing_ns(self, append_ns, contig, warn=True, replace_char='_'):
		""" Replaces "rightmost overhanging Ns" in trailing repeat with - character
		    No need to replace Ns in spacer: would mean that entire repeat was Ns and these elements are dropped
		"""
		repeat_end = self.items[-1].repeat_start + self.items[-1].repeat_length - 1
		contig_end = contig.length - append_ns
		self.items[-1].right_overhang = max(repeat_end - contig_end, 0)
		if self.items[-1].right_overhang > 0:
			if warn: 
				sys.stderr.write("Warning: Trailing Ns for array\n")
				sys.stderr.write("  before:%s\n" % self.items[-1].repeat_seq)
			seq = list(self.items[-1].repeat_seq)
			for i in range(1, self.items[-1].right_overhang+1):
				seq[-i] = replace_char
			self.items[-1].repeat_seq = ''.join(seq)
			if warn: 
				sys.stderr.write("   after:%s\n\n" % self.items[-1].repeat_seq)
			
	def fix_coords(self, append_ns):
		""" Because Ns are added to contigs, coorinates don't match up with genome
			Need to account for Ns added to contig and Ns leftover in repeat
		"""
		for item in self.items:
			item.repeat_start = item.repeat_start - append_ns + item.left_overhang
			
	def is_truncated(self, append_ns, contig, min_dist=50):
		""" Determine if array is likely truncated at the left or right side 
		"""
		self.start_pos = self.items[0].repeat_start
		self.end_pos = self.items[-1].repeat_start + self.items[-1].repeat_length - 1
		if self.start_pos - 1 <= min_dist: left = 'Yes'
		else: left = 'No'
		if len(contig.seq) - append_ns - self.end_pos - 1 <= min_dist: right = 'Yes'
		else: right = 'No'
		if right == 'Yes' and left == 'Yes':
			self.truncated = 'both'
		elif left == 'Yes':
			self.truncated = 'left'
		elif right == 'Yes':
			self.truncated = 'right'
		else:
			self.truncated = 'neither'
			
	def define_consensus(self):
		""" Define consensus direct repeat
		"""
		self.consensus_repeat = ""
		percent_id = []
		for i in range(len(self.items[0].repeat_seq)):
			all_chars = [item.repeat_seq[i] for item in self.items]
			count_gaps = all_chars.count('-')
			count_non_missing = len(all_chars) - all_chars.count('N')
			if float(count_gaps)/count_non_missing > 0.50:
				continue
			else:
				good_chars = [char for char in all_chars if char not in ['-', 'N']]
				good_count = len(good_chars)
				cons_base, cons_count = Counter(good_chars).most_common(1)[0]
				self.consensus_repeat += cons_base
				percent_id.append(100 * cons_count / float(good_count))
		self.percent_id = round(np.mean(percent_id), 2)

class ArrayElement:
	def __init__(self):
		self.repeat_seq = None
		self.spacer_seq = None
		self.repeat_start = None
		self.left_overhang = 0
		self.right_overhang = 0
	
class Contig:
	def __init__(self):
		id = None
		start = None
		end = None
		seq = None
		
def parse_args():
	import argparse	
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		add_help=False
		)
	parser.add_argument("-i", dest="inpath", metavar="FASTA", required=True,
						help="input FASTA file")
	parser.add_argument("-o", dest="outdir", metavar="PATH",
						help="output directory", required=True)
	parser.add_argument("-n", dest="append_ns", metavar="INT", type=int, default=50,
						help="number of Ns to append to beggining and end of contigs (50)")
	parser.add_argument("-b", dest="max_bp", metavar="INT", type=int, default=int(5e7), 
						help="max number of base-pairs to search at once (50,000,000)")
	parser.add_argument("-m", dest="max_seqs", metavar="INT", type=int, default=float('Inf'), 
						help="max number of sequences to search from input (search all)")
	parser.add_argument("-x", dest="xmx", metavar="INT", type=int, default=2000,
						help="Max memory in MB for Java (2000)")						
	parser.add_argument("-h", "--help", action="help", 
						help="show this help message and exit")		
	args = vars(parser.parse_args())

	return args

def fasta_offsets(fasta, max_bp, max_seqs):
	"""
	Split multi-fasta file into chunks, each with up to maxbp
	Return the offset of split start and the number of sequences in the split
	Useful to prevent CRISPR programs from crashing when fasta contains millions of contigs
	"""
	splits = []
	curbp = 0
	last_offset = 0
	cur_offset = 0 
	nseqs = 0
	tseqs = 0
	for line in open(fasta):
		if line[0] != '>':
			curbp += len(line.rstrip())
		else:
			nseqs += 1
			tseqs += 1
		if curbp >= max_bp:
			splits.append([last_offset, nseqs])
			nseqs = 0
			curbp = 0
			last_offset = cur_offset 	
		cur_offset += len(line)
		if tseqs >= max_seqs:
			return splits
	if nseqs > 0:
		splits.append([last_offset, nseqs])
	return splits

def write_tmp_seqs(fasta, outdir, offset, nseqs, append_ns, concat):
	"""
	Open fasta, seek to offset, read nseqs from fasta, append Ns to start & end, optionally concat seqs.
	Write sequences with unique ids to temporary directory (outdir/temp/seqs).
	Finally, return path to sequences, and a list of their ids, start/stop positions (important if concatenating),
	and their sequences.
	"""

	from Bio.SeqIO import parse 
	import uuid
	contigs = []
	seqs_file = open(fasta)
	seqs_file.seek(offset)
	seqs_bio = parse(seqs_file, 'fasta')
	
	tmp_path = '%s/temp/seqs/%s' % (outdir, str(uuid.uuid4()))
	with open(tmp_path, 'w') as tmp_file:
		if concat:
			tmp_file.write('>tmpfile\n')
		for i in range(nseqs):
			s = next(seqs_bio)
			contig = Contig()
			contig.id = s.id
			contig.seq = ('N'*append_ns) + str(s.seq).replace('X','N') + ('N'*append_ns)
			contig.length = len(contig.seq)
			if not concat:
				tmp_file.write('>%s\n' % contig.id)
				tmp_file.write('%s\n' % contig.seq)
				contig.start = 1
				contig.end = len(contig.seq)
			else:
				tmp_file.write(contig.seq)
				contig.start = 1 if len(contigs) == 0 else contigs[-1].end + 1
				contig.end = contig.start + len(contig.seq) - 1
			contigs.append(contig)
		if concat:
			tmp_file.write('\n')
			
	return tmp_path, contigs
	
def run_pipeline(args, offsets, append_ns, program):
	arrays = []
	
	for offset, nseqs in offsets:
		
		# write temp seqs
		concat_seqs = True if program == 'crt' else False
		tmp_path, contigs = write_tmp_seqs(args['inpath'], args['outdir'], offset, nseqs, append_ns, concat=concat_seqs)
		
		# run programs
		if program == 'crt':
			out, err, code = run_crt(tmp_path, args['tempdir'], xmx=args['xmx'])
		elif program == 'pilercr':
			out, err, code = run_pilercr(tmp_path, args['tempdir'])
			
		# check exit code
		if code != 0:
			sys.stderr.write("Warning: skipping sequences due to error: %s\n" % err)
			continue
			
		# parse output
		if program == 'crt':	
			my_arrays = parse_crt(out, contigs)
		elif program == 'pilercr':
			my_arrays = parse_pilercr(out, contigs)

		# clean up arrays; store them; remove temp files
		clean_up_arrays(my_arrays, append_ns, contigs)
		arrays += my_arrays

	return arrays

def clean_up_arrays(arrays, append_ns, contigs):
	contigs = dict([(_.id, _) for _ in contigs])
	for array in arrays:
		array.prune(max_nfraction=0.80, warn=False)
		array.replace_leading_ns(append_ns, contigs[array.contig_id], warn=False)
		array.replace_trailing_ns(append_ns, contigs[array.contig_id], warn=False)
		array.fix_coords(append_ns)
		array.is_truncated(append_ns, contigs[array.contig_id])

def write_crispr_files(arrays, outbase):
	
	last_contig = None
	array_num = 0
	
	array_file = open(outbase+'.arrays', 'w')
	array_header = ['contig_id', 'array_num', 'start_pos', 'end_pos', 'truncated',
		            'count_spacers', 'percent_id', 'repeat_length', 'consensus_repeat']
	array_file.write('\t'.join(array_header)+'\n')

	spacers_file = open(outbase+'.spacers', 'w')
	spacers_header = ['contig_id', 'array_num', 'repeat_num', 'start_pos', 'repeat_seq', 'spacer_seq']
	spacers_file.write('\t'.join(spacers_header)+'\n')

	for array in arrays:
	
		if array.contig_id != last_contig:
			array_num = 1
			last_contig = array.contig_id
		else:
			array_num += 1
		
		array_record = [
			array.contig_id,
			array_num, 
			array.start_pos,
			array.end_pos,
			array.truncated,
			len(array.items) - 1,
			array.percent_id,
			len(array.consensus_repeat),
			array.consensus_repeat
		]
		
		array_file.write('\t'.join([str(_) for _ in array_record])+'\n')
		
		for repeat_num, item in enumerate(array.items):
			spacers_record = [
				array.contig_id,
				array_num, 
				repeat_num+1, 
				item.repeat_start, 
				item.repeat_seq, 
				item.spacer_seq
			]
			spacers_file.write('\t'.join([str(_) for _ in spacers_record])+'\n')
			
	array_file.close()
	spacers_file.close()
									
if __name__ == "__main__":
	
	args = parse_args()
	
	args['tempdir'] = args['outdir']+'/temp'
	if not os.path.isdir(args['tempdir']):
		os.makedirs(args['tempdir']+'/seqs')
	
	offsets = fasta_offsets(args['inpath'], args['max_bp'], args['max_seqs'])
	
	for program in ['crt', 'pilercr']:
		arrays = run_pipeline(args, offsets, append_ns=50, program=program)
		write_crispr_files(arrays, outbase=args['outdir']+'/'+program)

	shutil.rmtree(args['tempdir'])

