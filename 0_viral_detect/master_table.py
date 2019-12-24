
import os, csv, sys

hmm_hits = sys.argv[1]
virfinder = sys.argv[2]
strand_switch = sys.argv[3]

info = {}

for r in csv.DictReader(open(hmm_hits), delimiter='\t'):
	info[r['contig_id']] = r
	info[r['contig_id']]['gene2fam'] = r['annotations']
	info[r['contig_id']]['vfr_pvalue'] = None
	info[r['contig_id']]['vfr_score'] = None

with open(virfinder) as f:
	next(f)
	for l in f:
		r = l.rstrip('\n').replace('"', '').split('\t')
		id = r[1].split()[0]
		info[id]['vfr_pvalue'] = r[-1]
		info[id]['vfr_score'] = r[-2]
		if 'multi' in r[1]:
			info[id]['depth'] = r[1].split()[2].split('=')[1]
		elif 'cov' in r[1]:
			info[id]['depth'] = r[1].split('_')[-1]
		else:
			info[id]['depth'] = None

for r in csv.DictReader(open(strand_switch), delimiter='\t'):
	info[r['contig_id']].update(r)
	info[r['contig_id']]['gene_len'] = r['avg_len']


fields = ['contig_id', 'length', 'genes', 'gene_len', 'cds_density', 'switch_rate', 'vfr_score', 'vfr_pvalue', 'vpfs', 'pfams', 'gene2fam']
sys.stdout.write('\t'.join(fields)+'\n')
for values in info.values():
	sys.stdout.write('\t'.join([str(values[f]) for f in fields])+'\n')

