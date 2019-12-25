
import Bio.SeqIO, numpy as np, sys, gzip

# initialize contigs
data = {}
file = gzip.open(sys.argv[1]) if sys.argv[1].split('.')[-1]=='gz' else open(sys.argv[1])
for r in Bio.SeqIO.parse(file, 'fasta'):
	contig_id = r.id
	data[contig_id] = {}
	data[contig_id]['contig_id'] = contig_id
	data[contig_id]['length'] = len(str(r.seq))
	data[contig_id]['strand'] = []
	data[contig_id]['genes'] = set([])
	data[contig_id]['gene_length'] = []
file.close()

# add gene counts
file = gzip.open(sys.argv[2]) if sys.argv[2].split('.')[-1]=='gz' else open(sys.argv[2])
for r in Bio.SeqIO.parse(file, 'fasta'):
	gene_id = r.id
	contig_id = r.id.rsplit('_', 1)[0]
	start = int(r.description.split()[2])
	stop = int(r.description.split()[4])
	length = (stop - start + 1)
	strand = r.description.split()[6]
	data[contig_id]['strand'].append(strand)
	data[contig_id]['genes'].add(gene_id)
	data[contig_id]['gene_length'].append(length)

# summary statistics
for id in data:
	num_genes = len(data[id]['genes'])
	switches = 0
	strands = data[id]['strand']
	for i in range(len(strands)-1):
		if strands[i] != strands[i+1]:
			switches += 1
	data[id]['num_genes'] = len(data[id]['genes'])
	data[id]['switch_rate'] = round(100.0*switches/num_genes,2) if num_genes > 0 else 0.0
	data[id]['avg_len'] = round(np.mean(data[id]['gene_length']),2) if num_genes > 0 else 0.0
	data[id]['cds_density'] = round(100.0*sum(data[id]['gene_length'])/data[id]['length'], 2)

# write results
fields = ['contig_id', 'length', 'num_genes', 'avg_len', 'switch_rate', 'cds_density']
sys.stdout.write('\t'.join(fields)+'\n')
for values in data.values():
	sys.stdout.write('\t'.join([str(values[f]) for f in fields])+'\n')













