#!/usr/bin/env python

import os, csv, operator, argparse

# fetch args
parser = argparse.ArgumentParser()
parser.add_argument('--in_features', required=True, metavar='PATH')
parser.add_argument('--out_features', required=True, metavar='PATH')
parser.add_argument('--in_seqs_base', required=True, metavar='PATH')
parser.add_argument('--out_seqs_base', required=True, metavar='PATH')
args = vars(parser.parse_args())

# read parameters
rule_sets = {}
for r in csv.DictReader(open('classification_rules.tsv'), delimiter='\t'):
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
lengths = rule_sets.keys()
rows = []
for r in csv.DictReader(open(args['in_features']), delimiter='\t'):
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

# write tabular output
with open(args['out_features'], 'w') as out:
	fields = ['contig_id', 'length',
			  'genes', 'gene_len', 'cds_density', 'switch_rate',
			  'vfr_score', 'vfr_pvalue',
			  'vpfs', 'vpf_pct', 'pfams', 'pfam_pct', 'gene2fam']
	out.write('\t'.join(fields)+'\n')
	for row in rows:
		row = [str(row[f]) for f in fields]
		out.write('\t'.join(row)+'\n')

# write sequences
import Bio.SeqIO
contigs = set([r['contig_id'] for r in rows])
for type in ['fna', 'faa', 'ffn']:
	inpath = args['in_seqs_base']+'.'+type
	outpath = args['out_seqs_base']+'.'+type
	out = open(outpath, 'w')
	for r in Bio.SeqIO.parse(inpath, 'fasta'):
		c = r.id if type == 'fna' else r.id.rsplit('_', 1)[0]
		if c in contigs:
			out.write('>'+r.description+'\n'+str(r.seq)+'\n')



