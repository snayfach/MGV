#!/bin/bash -l

import os, Bio.SeqIO, subprocess as sp, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--out_dir', required=True, metavar='PATH')
parser.add_argument('--threads', type=int, default=1, metavar='INT')
args = vars(parser.parse_args())

# make dir
if not os.path.exists(args['out_dir']):
	os.makedirs(args['out_dir'])

# read proteins
faa = dict([(r.id, str(r.seq).rstrip('*')) for r in Bio.SeqIO.parse(args['in_faa'], 'fasta')])

print ("make diamond db")
m8_path = '%s/dmd.tsv' % args['out_dir']
cmd = "diamond makedb --in %s --db %s/db " % (args['in_faa'], args['out_dir'])
cmd += "--threads %s" % args['threads']
p = sp.Popen(cmd, shell=True)
p.wait()

print ("run diamond")
cmd = "diamond blastp --query %s " % args['in_faa']
cmd += "--db %s/db " % args['out_dir']
cmd += "--out %s/dmnd.tsv " % args['out_dir']
cmd += "--threads %s " % args['threads']
cmd += "--max-target-seqs 10000 --evalue 1e-5 "
cmd += "--query-cover 50 --subject-cover 50"
p = sp.Popen(cmd, shell=True)
p.wait()

print ("run mcl")
mcl_input = '%s/mcl_input.tsv' % args['out_dir']
mcl_output = '%s/mcl_output.tsv' % args['out_dir']
with open(mcl_input, 'w') as f:
	for line in open('%s/dmnd.tsv' % args['out_dir']):
		r = line.split()
		if r[0] != r[1]:
			row = [r[0], r[1], r[-1]]
			f.write('\t'.join(row)+'\n')
cmd = "mcl %s --abc -I 1.4 -te %s -o %s" % (mcl_input, args['threads'], mcl_output)
p = sp.Popen(cmd, shell=True)
p.wait()

print ("select protein clusters (PCs)")
pcs = {}
genome_ids = set([_[1:].split()[0].rsplit('_', 1)[0] for _ in open(args['in_faa']) if _[0]=='>'])
for pc_id, line in enumerate(open(mcl_output)):
	pc_id = str(pc_id)
	pc_genes = line.split()
	pc_otus = list(set([_.split('_')[0] for _ in pc_genes]))
	prev = 100.0*len(pc_otus)/len(genome_ids)
	copy = 1.0*len(pc_genes)/len(pc_otus)
	pcs[pc_id] = {'genes':pc_genes, 'prev':prev, 'copy':copy}
marker_ids = []
for id, values in pcs.items():
	if values['prev'] > 10 and values['copy'] < 1.1:
		marker_ids.append(id)
open(args['out_dir']+'/marker_ids.list', 'w').write('\n'.join(marker_ids)+'\n')

print ("write seqs for each PC")
faa_dir = '%s/faa' % args['out_dir']
if not os.path.exists(faa_dir):
	os.makedirs(faa_dir)
for id in marker_ids:
	out = open('%s/%s.faa' % (faa_dir, id), 'w')
	for gene_id in pcs[id]['genes']:
		out.write('>'+gene_id+'\n'+faa[gene_id]+'\n')

print ("build multiple seq alignment (MSA) for each PC, round 1")
msa_dir = '%s/msa' % args['out_dir']
if not os.path.exists(msa_dir):
	os.makedirs(msa_dir)
for file in os.listdir(faa_dir):
	inpath = faa_dir+'/'+file
	outpath = msa_dir+'/'+file
    cmd = "famsa -t %s %s %s" % (args['threads'], inpath, outpath)
	p = sp.Popen(cmd, shell=True)
	p.wait()

print ("build hmm from each MSA")
hmm_dir = '%s/hmm' % args['out_dir']
if not os.path.exists(hmm_dir):
	os.makedirs(hmm_dir)
for file in os.listdir(msa_dir):
	inpath = msa_dir+'/'+file
	outpath = hmm_dir+'/'+file.replace('.faa', '.hmm')
	cmd = "hmmbuild %s %s" % (outpath, inpath)
	p = sp.Popen(cmd, shell=True)
	p.wait()
all_hmms = '%s/all.hmms' % args['out_dir']
with open(all_hmms, 'w') as f:
	for file in os.listdir(hmm_dir):
		for line in open(hmm_dir+'/'+file):
			f.write(line)

print ("search hmms versus original proteins")
cmd = "hmmsearch --noali --cpu %s " % args['threads']
cmd += "--tblout %s/hmmsearch.txt " % args['out_dir']
cmd += "%s/all.hmms " % args['out_dir']
cmd += "%s > /dev/null" % args['in_faa']
p = sp.Popen(cmd, shell=True)
p.wait()

print ("extract proteins for top hmm hits")
hits_dir = '%s/hits' % args['out_dir']
if not os.path.exists(hits_dir):
	os.makedirs(hits_dir)
hits = {}
for line in open('%s/hmmsearch.txt' % args['out_dir']):
	if line[0]=='#': continue
	r = line.split()
	if float(r[4]) > 1e-5: continue
	elif r[0] not in hits:
		hits[r[0]] = r
	elif float(r[5]) > float(hits[r[0]][5]):
		hits[r[0]] = r
marker_to_hits = {}
for r in hits.values():
	marker_id = r[2]
	if marker_id not in marker_to_hits:
		marker_to_hits[marker_id] = []
	gene = faa[r[0]]
	marker_to_hits[marker_id].append([r[0], gene])
for marker_id in marker_to_hits:
	outpath = '%s/%s.faa' % (hits_dir, marker_id)
	out = open(outpath, 'w')
	for gene_id, seq in marker_to_hits[marker_id]:
		out.write('>'+gene_id+'\n'+seq+'\n')

print ("build multiple seq alignment (MSA) for each PC, round 2")
msa2_dir = '%s/msa2' % args['out_dir']
if not os.path.exists(msa2_dir):
	os.makedirs(msa2_dir)
for file in os.listdir(hits_dir):
	inpath = hits_dir+'/'+file
	outpath = msa2_dir+'/'+file
	cmd = "famsa -t %s %s %s" % (args['threads'], inpath, outpath)
	p = sp.Popen(cmd, shell=True)
	p.wait()

print ("trim multiple seq alignments")
trim2_dir = '%s/trim2' % args['out_dir']
if not os.path.exists(trim2_dir):
	os.makedirs(trim2_dir)
for file in os.listdir(msa2_dir):
	inpath = msa2_dir+'/'+file
	outpath = trim2_dir+'/'+file
	cmd = "trimal -gt 0.5 -in %s -out %s" % (inpath, outpath)
	p = sp.Popen(cmd, shell=True)
	p.wait()

print ("concatenate alignments and fill missing genes with gaps")
msa_catted = dict([(id, '') for id in genome_ids])
for index, file in enumerate(os.listdir(trim2_dir)):
	marker_id = file.split('.')[0]
	genome_to_gene = {}
	inpath = '%s/%s' % (trim2_dir, file)
	for r in Bio.SeqIO.parse(inpath, format='fasta'):
		genome_id = r.id.rsplit('_', 1)[0]
		genome_to_gene[genome_id] = r
	marker_length = len(list(genome_to_gene.values())[0].seq)
	for genome_id in genome_ids:
		if genome_id not in genome_to_gene:
			msa_catted[genome_id] += "-" * marker_length
		else:
			gene = genome_to_gene[genome_id]
			msa_catted[genome_id] += str(gene.seq).replace('X', '-')

print ("count gaps per genome")
genome_to_gaps = {}
for genome_id in genome_ids:
	seq = msa_catted[genome_id]
	perc_gaps = round(100*seq.count('-')/float(len(seq)), 2)
	genome_to_gaps[genome_id] = perc_gaps

print ("write concatenated alignment for genomes with <90% gaps")
with open('%s/concat.faa' % args['out_dir'], 'w') as f:
	for genome_id, seq in msa_catted.items():
		perc_gaps = genome_to_gaps[genome_id]
		if perc_gaps < 90:
			f.write('>%s percent_gaps=%s\n' % (genome_id, perc_gaps))
			f.write('%s\n' % seq)

print ("built phylogeny with FastTree")
my_env = os.environ
my_env["OMP_NUM_THREADS"] = "16"
cmd = "FastTreeMP %s/concat.faa > %s/concat.tree" % (args['out_dir'], args['out_dir'])
p = sp.Popen(cmd, shell=True, env=my_env)
p.wait()
