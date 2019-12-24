
import os

parser = argparse.ArgumentParser()
parser.add_argument('--in_faa', required=True, metavar='PATH')
parser.add_argument('--out_base', required=True, metavar='PATH')
parser.add_argument('--threads', default=1, type=int, metavar='INT')
parser.add_argument('--min_gene_sharing', default=50, type=float, metavar='FLOAT')
args = vars(parser.parse_args())


query_lengths = {}
with open('contig_length.tsv') as file:
	next(file)
	for line in file:
		query, length = line.split()[0:2]
		query_lengths[query] = int(length)

target_lengths = {}
with open('../0_input_data/input_lengths.tsv') as file:
	next(file)
	for line in file:
		target, length = line.split()[0:2]
		target_lengths[target] = int(length)

rows = []

indir = '../3_aaid/aaid/'
for index, file in enumerate(os.listdir(indir)):
	print index
	with open(indir+file) as handle:
		next(handle)
		for line in handle:
			query, target, num_genes, hits, cov, aaid = line.split()
			if float(cov) < 50.0 or float(aaid) < 80.0:
				continue
			qlength = query_lengths[query]
			tlength = target_lengths[target]
			complete = min([100.0*qlength/tlength, 99.0])
			row = [query, target, num_genes, hits, cov, aaid, qlength, tlength, complete]
			rows.append(row)

out = open('completeness.tsv', 'w')
header = ['query', 'target', 'num_genes', 'hits', 'cov', 'aaid', 'qlength', 'tlength', 'complete']
out.write('\t'.join(header)+'\n')
for row in rows:
	out.write('\t'.join([str(_) for _ in row])+'\n')




