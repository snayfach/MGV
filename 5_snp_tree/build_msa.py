#!/usr/bin/env python

import subprocess as sp, os, time, sys
import numpy as np

def fetch_args():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--in', type=str, required=True,
		help="""Input directory; same as --out from align_genomes.py""")
	parser.add_argument('--out', type=str, required=True,
		help="""Output path""")
	parser.add_argument('--max_genomes', type=int, default=None,
		help="""Maximum # of genomes to process; useful for quick tests""")
	parser.add_argument('--max_sites', type=int, default=None,
		help="""Maximum # of sites to process; useful for quick tests""")
	parser.add_argument('--max_gaps_col', type=float, default=100,
		help="""Max fraction gaps per position (100)""")
	parser.add_argument('--max_gaps_seq', type=float, default=100,
		help="""Max fraction of gaps per genome (100)""")
	return vars(parser.parse_args())

def parse_seqs(path):
	with open(path) as file:
		try: id = next(file).split()[0].lstrip('>')
		except: return
		seq = ''
		for line in file:
			if line[0]=='>':
				yield id, seq
				try: id = line.split()[0].lstrip('>')
				except: return
				seq = ''
			else:
				seq += line.rstrip()
		yield id, seq

def parse_coords(fpath):
	fields = [('s1',int),('e1',int),
			  ('s2',int),('e2',int),
	          ('len1',int),('len2',int),
	          ('pid',float),
	          ('c1',str),('c2',str)]
	with open(fpath) as f:
		for i in range(5):
			next(f)
		for l in f:
			values = l.replace(' | ', ' ').split()
			yield dict([(f[0],f[1](v)) for f,v in zip(fields, values)])

def parse_snps(fpath):
	fields = [('p1',int),('b1',str),('b2',str),('p2',int),
			  ('buf',int),('dist',int),
	          ('r',int),('q',int),
	          ('s1',int),('s2',int),
	          ('c1',str),('c2',str)]
	with open(fpath) as f:
		for i in range(5):
			next(f)
		for l in f:
			values = l.replace(' | ', ' ').split()
			yield dict([(f[0],f[1](v)) for f,v in zip(fields, values)])

def keep_col(x, max_gaps_col):
	if 100.0*(x == '-').sum()/float(len(x)) >= max_gaps_col:
		return False
	else:
		return True

def is_snp(x):
	y = x[x != '-']
	if (y != y[0]).sum() > 0:
		return True
	else:
		return False

if __name__ == "__main__":

	start = time.time()
	
	args = fetch_args()
	args['aln_dir'] = os.path.join(args['in'], 'aln')
	
	if not os.path.exists(args['in']):
		sys.exit("Error: dir does not exist: %s" % args['in'])
	
	print("Reading reference genome")
	ref = {}
	for id, seq in parse_seqs(os.path.join(args['in'], 'reference.fna')):
		ref[id] = np.array(list(seq.upper()))
	print("   count contigs: %s" % len(ref))
	print("   count sites: %s" % sum([len(_) for _ in ref.values()]))

	print("Initializing alignments")
	genome_ids = os.listdir(args['aln_dir'])
	if args['max_genomes']:
		genome_ids = genome_ids[:args['max_genomes']]
	print("   count genomes: %s" % len(genome_ids))
	genomes = {}
	for index, genome_id in enumerate(genome_ids):
		genomes[genome_id] = {}
		for id, seq in ref.items():
			genomes[genome_id][id] = np.array(['-']*len(seq))
	# read alignments
	for index, genome_id in enumerate(genome_ids):
		fpath = '%s/%s/coords' % (args['aln_dir'], genome_id)
		aln_length = 0
		for r in parse_coords(fpath):
			aln_length += (r['e1'] - r['s1'])
			genomes[genome_id][r['c1']][r['s1']-1:r['e1']] = ref[r['c1']][r['s1']-1:r['e1']]
	# snps
	for index, genome_id in enumerate(genome_ids):
		fpath = '%s/%s/snps' % (args['aln_dir'], genome_id)
		for r in parse_snps(fpath):
			if r['b1'] == '.':
				continue
			elif r['b2'] == '.':
				genomes[genome_id][r['c1']][r['p1']-1] = '-'
			else:
				genomes[genome_id][r['c1']][r['p1']-1] = r['b2']
	# concatenate contigs
	for index, genome_id in enumerate(genomes):
		genomes[genome_id] = np.concatenate([seq for seq in genomes[genome_id].values()])
	if args['max_sites']:
		for genome_id in genomes:
			genomes[genome_id] = genomes[genome_id][:args['max_sites']]

	print("Pruning low frequency sites")
	print("   max gaps per col: %s%%" % args['max_gaps_col'])
	tmpdata = np.array([seq for seq in genomes.values()])
	sites = np.apply_along_axis(keep_col, 0, tmpdata, args['max_gaps_col'])
	for genome_id in genomes:
		genomes[genome_id] = genomes[genome_id][sites]
	print("   retained sites: %s" % sites.sum())

	print("Pruning invariant sites")
	tmpdata = np.array([seq for seq in genomes.values()])
	snps = np.apply_along_axis(is_snp, 0, tmpdata)
	for genome_id in genomes:
		genomes[genome_id] = genomes[genome_id][snps]
	print("   retained sites: %s" % snps.sum())

	print("Pruning genomes with too many gaps")
	print("   max gaps per seq: %s%%" % args['max_gaps_seq'])
	remove = []
	for genome_id, seq in genomes.items():
		fract_missing = 100.0*(seq == '-').sum()/float(len(seq))
		if fract_missing > args['max_gaps_seq']:
			remove.append(genome_id)
	for genome_id in remove:
		del genomes[genome_id]
	print("   retained genomes: %s" % len(genomes))

	print("Writing fasta")
	print("   path: %s" % args['out'])
	with open(args['out'], 'w') as f:
		for genome_id, seq in genomes.items():
			percent_missing = 100*(seq == '-').sum()/float(len(seq))
			f.write('>%s percent_missing=%s\n' % (genome_id, round(percent_missing,2)))
			f.write(''.join(seq)+'\n')







