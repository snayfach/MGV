
import subprocess as sp, os, numpy as np, sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--out_base', required=True, metavar='PATH')
parser.add_argument('--threads', default=1, type=int, metavar='INT')
parser.add_argument('--min_gene_sharing', default=50, type=float, metavar='FLOAT')
args = vars(parser.parse_args())

print("perform search...")
cmd = "diamond blastp "
cmd += "--threads %s " % args['threads']
cmd += "--db input/viral_proteins "
cmd += "--query %s " % args['in_faa']
cmd += "--out %s.blastout " % args['out_base']
cmd += "--outfmt 6 "
cmd += "--evalue 1e-5 "
cmd += "--max-target-seqs 10000 "
p = sp.Popen(cmd, shell=True)
p.wait()

print("init seqs...")
from collections import defaultdict
data = {}
for line in open(args['in_faa']):
	if line[0]=='>':
		contig_id = line[1:].split()[0].rsplit('_', 1)[0]
		if contig_id not in data:
			data[contig_id] = {'genes':0, 'alns':defaultdict(list)}
		data[contig_id]['genes'] += 1

print("store alignments...")
gene_to_genome = dict([line.split()[0:2] for line in open("input/gene_to_genome.tsv")])
for line in open("%s.blastout" % args['out_base']):
	aln = line.split()
	qgene, tgene = aln[0:2]
	qcontig = qgene.rsplit('_', 1)[0]
	tcontig = gene_to_genome[tgene]
	data[qcontig]['alns'][tcontig].append(aln)

print("compute aaid...")
rows = []
for query in data:
	num_genes = data[query]['genes']
	bts = {}
	for target in data[query]['alns']:
		# store best hits to target genes
		bhs = {}
		for aln in data[query]['alns'][target]:
			qgene, tgene = aln[0:2]
			pid = float(aln[2])
			if qgene not in bhs:
				bhs[qgene] = (tgene, pid)
			elif pid > bhs[qgene][1]:
				bhs[qgene] = (tgene, pid)
		# compute aaid between query and target
		hits = len(bhs)
		cov = 100.0*hits/num_genes
		aaid = np.mean([_[1] for _ in bhs.values()])
		# store aaid to top target
		if cov < args['min_gene_sharing']:
			continue
		if target not in bts:
			bts[target] = (aaid, cov, hits)
		elif aaid > bts[target][0]:
			bts[target] = (aaid, cov, hits)

	if len(bts) > 0:
		target = list(bts.keys())[0]
		aaid, cov, hits = bts[target]
		rows.append([query, target, num_genes, hits, cov, aaid])

print("write results...")
out = open("%s.aai" % args['out_base'], "w")
header = ['query', 'target', 'num_genes', 'hits', 'cov', 'aaid']
out.write('\t'.join(header)+'\n')
for row in rows:
	out.write('\t'.join([str(_) for _ in row])+'\n')

#print "clean up"
#os.remove("search/%s.m8" % file_num)









