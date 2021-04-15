#!/usr/bin/env python

import os, csv, operator, argparse

# fetch args
parser = argparse.ArgumentParser()
parser.add_argument('--features', required=True, metavar='PATH')
parser.add_argument('--in_base', required=True, metavar='PATH')
parser.add_argument('--out_base', required=True, metavar='PATH')
args = vars(parser.parse_args())

# read parameters
rule_sets = {}
for r in csv.DictReader(open('input/classification_rules.tsv'), delimiter='\t'):
	length = int(r['length'])
	if length not in rule_sets:
		rule_sets[length] = []
	rule_sets[length].append([
		('genes', operator.ge, int(r['min_genes'])),
		('vpf_pct', operator.ge, float(r['min_vpf_pct'])),
		('pfam_pct', operator.le, float(r['max_pfam_pct'])),
		('vfr_pvalue', operator.le, float(r['max_pvalue'])),
		('switch_rate', operator.le, float(r['max_switch']))])

# read data
lengths = list(rule_sets.keys())
rows = []
for r in csv.DictReader(open(args['features']), delimiter='\t'):
	# filter row
	r['length'] = int(r['length'])
	if r['length'] < 1000:
		continue
	# format row
	r['genes'] = int(r['genes'])
	r['vpf_pct'] = 100.0*int(r['vpfs'])/int(r['genes']) if int(r['genes']) > 0 else 0.0
	r['pfam_pct'] = 100.0*int(r['pfams'])/int(r['genes']) if int(r['genes']) > 0 else 0.0
	r['vfr_pvalue'] = float(r['vfr_pvalue']) if r['vfr_pvalue'] != 'None' else 1.0
	r['switch_rate'] = float(r['switch_rate'])
	# determine rule set to use
	deltas = [abs(l - r['length']) for l in lengths]
	rule_set = rule_sets[lengths[deltas.index(min(deltas))]]
	# apply rule set
	bools = []
	for rule in rule_set:
		bools.append(all([my_operator(r[field], value) for field, my_operator, value in rule]))
	if any(bools):
		rows.append(r)

# viral contigs
viral_contigs = set([r['contig_id'] for r in rows])

# write tabular output
with open(args['out_base']+'.tsv', 'w') as out:
	fields = ['contig_id', 'length',
			  'genes', 'gene_len', 'cds_density', 'switch_rate',
			  'vfr_score', 'vfr_pvalue',
			  'vpfs', 'vpf_pct', 'pfams', 'pfam_pct', 'gene2fam']
	out.write('\t'.join(fields)+'\n')
	for row in rows:
		row = [str(row[f]) for f in fields]
		out.write('\t'.join(row)+'\n')

# write gff
inpath = args['in_base']+'.gff'
if os.path.exists(inpath):
	outpath = args['out_base']+'.gff'
	out = open(outpath, 'w')
	handle = open(inpath)
	out.write(next(handle))
	write_flag = False
	for line in handle:
		if '# Sequence Data:' in line:
			id = line.split('seqhdr=')[1].split()[0].lstrip('"')
			if id in viral_contigs: write_flag = True
			else: write_flag = False
		if write_flag: out.write(line)

# write sequences
import Bio.SeqIO
for type in ['fna', 'faa', 'ffn']:
	inpath = args['in_base']+'.'+type
	outpath = args['out_base']+'.'+type
	out = open(outpath, 'w')
	for r in Bio.SeqIO.parse(inpath, 'fasta'):
		c = r.id if type == 'fna' else r.id.rsplit('_', 1)[0]
		if c in viral_contigs:
			out.write('>'+r.description+'\n'+str(r.seq)+'\n')



