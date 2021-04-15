#!/usr/bin/env python

import sys, os, gzip, Bio.SeqIO

MAX_EVALUE = 1e-10

# input arguments
fna_path = sys.argv[1]
faa_path = sys.argv[2]
vpf_path = sys.argv[3]
pfam_path = sys.argv[4]

# black-listed families
exclude = set([])
for file in ['microbial_vpfs.list', 'plasmid_vpfs.list', 'viral_pfams.list']:
	inpath = 'input/exclude_hmms/'+file
	exclude.update([_.rstrip() for _ in open(inpath)])

# init sequences
contigs = {}
file = gzip.open(fna_path) if fna_path.split('.')[-1] == 'gz' else open(fna_path)
for r in Bio.SeqIO.parse(file, 'fasta'):
	contigs[r.id] = {}
	contigs[r.id]['vpfs'] = []
	contigs[r.id]['pfams'] = []
	contigs[r.id]['length'] = len(str(r.seq))
	contigs[r.id]['genes'] = 0

# add gene counts
gene_length = {}
file = gzip.open(faa_path) if faa_path.split('.')[-1] == 'gz' else open(faa_path)
for r in Bio.SeqIO.parse(file, 'fasta'):
	contig_id = r.id.rsplit('_', 1)[0]
	contigs[contig_id]['genes'] += 1
	gene_length[r.id] = len(str(r.seq))
total_genes = len(gene_length)

# identify best hits by score
hits = {}
for dbname, hmmout in [['vpfs', vpf_path], ['pfams', pfam_path]]:
	file = gzip.open(hmmout) if hmmout.split('.')[-1] == 'gz' else open(hmmout)
	for l in open(hmmout):
		if l[0] == '#': continue
		elif len(l.split()) < 6: continue
		query = l.split()[0]
		score = float(l.split()[5])
		evalue = float(l.split()[4])#/total_genes
		if dbname == 'vpfs': target = l.split()[2]
		else: target = l.split()[3]
		if query not in hits:
			hits[query] = {'target':target, 'evalue':evalue, 'score':score, 'rec':l, 'db':dbname}
		elif score > hits[query]['score']:
			hits[query] = {'target':target, 'evalue':evalue, 'score':score, 'rec':l, 'db':dbname}

# count hits
for gene_id in hits:
	if hits[gene_id]['evalue'] > MAX_EVALUE:
		continue
	elif hits[gene_id]['target'] in exclude:
		continue

	contig_id = gene_id.rsplit('_', 1)[0]
	gene_num = gene_id.split('_')[-1]
	target_id = hits[gene_id]['target']

	if hits[gene_id]['db'] == 'vpfs':
		contigs[contig_id]['vpfs'].append(gene_num+':'+target_id)
	elif hits[gene_id]['db'] == 'pfams':
		contigs[contig_id]['pfams'].append(gene_num+':'+target_id)

# report results
header = ['contig_id', 'length', 'genes', 'vpfs', 'pfams', 'annotations']
sys.stdout.write('\t'.join(header)+'\n')
for contig_id in contigs:
	info = contigs[contig_id]
	annotations = ','.join(info['vpfs'] + info['pfams'])
	row = [contig_id, info['length'], info['genes'], len(info['vpfs']), len(info['pfams']), annotations]
	sys.stdout.write('\t'.join([str(_) for _ in row])+'\n')


