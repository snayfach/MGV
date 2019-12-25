#!/usr/bin/env python

import subprocess as sp, os, numpy as np, sys, Bio.SeqIO, csv, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--in_blast', required=True, metavar='PATH')
parser.add_argument('--out_tsv', required=True, metavar='PATH')
args = vars(parser.parse_args())

from collections import defaultdict
genome_size = defaultdict(int)
gene_to_genome = {}
for index, r in enumerate(Bio.SeqIO.parse(args['in_faa'], 'fasta')):
	genome_id = r.id.rsplit('_', 1)[0]
	genome_size[genome_id] += 1
	gene_to_genome[r.id] = genome_id
genome_alns = {}
for genome_id in genome_size:
	genome_alns[genome_id] = {}

print "parse"
for index, line in enumerate(open(args['in_blast'])):
	
	aln = line.split()
	gene = aln[0]
	query = gene_to_genome[aln[0]]
	target = gene_to_genome[aln[1]]
	score = float(aln[-1])
	
	if target not in genome_alns[query]:
		genome_alns[query][target] = {}
	if gene not in genome_alns[query][target]:
		genome_alns[query][target][gene] = aln
	elif score > float(genome_alns[query][target][gene][-1]):
		genome_alns[query][target][gene] = aln

print "compute"
rows = []
for query in genome_alns:
	query_genes = genome_size[query]
	for target in genome_alns[query]:
		target_genes = genome_size[target]
		alns = genome_alns[query][target].values()
		shared_genes = len(alns)
		qcov = 100.0*shared_genes/query_genes
		tcov = 100.0*shared_genes/target_genes
		aai = np.mean([float(_[2]) for _ in alns])
		row = [query, target, query_genes, target_genes, shared_genes, qcov, tcov, aai]
		rows.append(row)

print "write"
with open(args['out_tsv'], 'w') as out:
	header = ['qname', 'tname', 'qgenes', 'tgenes', 'sgenes', 'qcov', 'tcov', 'aai']
	out.write('\t'.join(header)+'\n')
	for row in rows:
		out.write('\t'.join([str(_) for _ in row])+'\n')



