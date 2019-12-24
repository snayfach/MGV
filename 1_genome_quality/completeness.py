
import subprocess as sp, os, numpy as np, sys, argparse
import Bio.SeqIO, csv
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--db_dir', required=True, metavar='PATH')
parser.add_argument('--in_fna', required=True, metavar='PATH')
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--out_m8', required=True, metavar='PATH')
parser.add_argument('--out_tsv', required=True, metavar='PATH')
parser.add_argument('--threads', default=1, type=int, metavar='INT')
parser.add_argument('--min_gene_sharing', default=50, type=float, metavar='FLOAT')
parser.add_argument('--min_aai', default=80, type=float, metavar='FLOAT')
args = vars(parser.parse_args())

print("perform search...")
cmd = "diamond blastp "
cmd += "--threads %s " % args['threads']
cmd += "--db %s/viral_proteins " % args['db_dir']
cmd += "--query %s " % args['in_faa']
cmd += "--out %s " % args['out_m8']
cmd += "--outfmt 6 "
cmd += "--evalue 1e-5 "
cmd += "--max-target-seqs 10000 "
p = sp.Popen(cmd, shell=True)
p.wait()

print("read query and target info...")
queries = {}
for r in Bio.SeqIO.parse(args['in_fna'], 'fasta'):
	queries[r.id] = {'query':r.id, 'dna':len(r.seq), 'genes':0, 'alns':defaultdict(list)}
	queries[r.id]['target'] = None
	queries[r.id]['aai'] = None
	queries[r.id]['shared_genes'] = None
	queries[r.id]['completeness'] = None
for r in Bio.SeqIO.parse(args['in_faa'], 'fasta'):
	queries[r.id.rsplit('_', 1)[0]]['genes'] += 1
targets = {}
with open('%s/genome_lengths.tsv' % args['db_dir']) as file:
	next(file)
	for line in file:
		target_id, length = line.split()[0:2]
		targets[target_id] = {'dna':int(length), 'genes':0}
gene_to_target = {}
for r in csv.DictReader(open('%s/gene_to_genome.tsv' % args['db_dir']), delimiter='\t'):
	gene_to_target[r['protein_id']] = r['contig_id']
	targets[r['contig_id']]['genes'] += 1

print("store alignments...")
for line in open(args['out_m8']):
	aln = line.split()
	qgene, tgene = aln[0:2]
	query = qgene.rsplit('_', 1)[0]
	target = gene_to_target[tgene]
	queries[query]['alns'][target].append(aln)

print("find nearest reference...")
for query in queries:
	# find best genome hit
	bts = ()
	for target in queries[query]['alns']:
		# store best hits to target genes
		bhs = {}
		for aln in queries[query]['alns'][target]:
			qgene, tgene = aln[0:2]
			pid = float(aln[2])
			if qgene not in bhs:
				bhs[qgene] = (tgene, pid)
			elif pid > bhs[qgene][1]:
				bhs[qgene] = (tgene, pid)
		# compute aai between query and target
		aai = np.mean([_[1] for _ in bhs.values()])
		shared = len(bhs)
		score = sum([_[1] for _ in bhs.values()])
		# store aai to genome
		if len(bts) == 0:
			bts = (target, aai, shared, score)
		elif score > bts[-1]:
			bts = (target, aai, shared, score)

	# store best hit
	if len(bts) > 0:
		target, aai, shared, score = bts
		queries[query]['target'] = target
		queries[query]['aai'] = aai
		queries[query]['shared_genes'] = shared

print("estimate completeness...")
rows = []
for query, rec in queries.items():
	# query, qlength, qgenes
	row = [query, rec['dna'], rec['genes']]
	# target, tlength, tgenes, tdatabase
	if rec['target'] is not None:
		if 'GCA' in rec['target']: db = 'NCBI'
		elif '____' in rec['target']: db = 'IMG/VR'
		else: db = 'UGV'
		row += [rec['target'], targets[rec['target']]['dna'], targets[rec['target']]['genes'], db]
	else:
		row += ['NA', 'NA', 'NA', 'NA']
	# shared_genes, aai, percent_shared, pass_shared, pass_aai, completeness
	if rec['target'] is not None:
		percent_shared = 100.0*rec['shared_genes']/rec['genes']
		pass_genes = 'Yes' if percent_shared >= args['min_gene_sharing'] else 'No'
		pass_aai = 'Yes' if rec['aai'] >= args['min_aai'] else 'No'
		complete = 'NA'
		if pass_genes=='Yes' and pass_aai=='Yes':
			complete = min([100.0*rec['dna']/targets[rec['target']]['dna'], 99.0])
		row += [rec['shared_genes'], rec['aai'], percent_shared, pass_genes, pass_aai, complete]
	else:
		row += ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
	rows.append(row)

out = open(args['out_tsv'], 'w')
fields = ['query', 'qlength', 'qgenes', 'target', 'tlength', 'tgenes', 'tdatabase', 'shared_genes', 'aai', 'percent_shared', 'pass_shared', 'pass_aai', 'completeness']
out.write('\t'.join(fields)+'\n')
for row in rows:
	out.write('\t'.join([str(_) for _ in row])+'\n')

