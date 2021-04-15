
import csv, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_aai', required=True, metavar='PATH')
parser.add_argument('--min_num_shared', type=int, required=True, metavar='INT')
parser.add_argument('--min_percent_shared', type=float, required=True, metavar='FLOAT')
parser.add_argument('--min_aai', type=float, required=True, metavar='FLOAT')
parser.add_argument('--out_tsv', required=True, metavar='PATH')
args = vars(parser.parse_args())

rows = []
for i, r in enumerate(csv.DictReader(open(args['in_aai']), delimiter='\t')):
	qcov, tcov = float(r['qcov']), float(r['tcov'])
	cov = min([qcov, tcov])
	aai = float(r['aai'])
	shared = int(r['sgenes'])
	score = cov/100.0 * aai/100.0
	if not (
			(cov>=args['min_percent_shared'] or shared>=args['min_num_shared'])
			and aai >= args['min_aai']):
		continue
	row = [r['qname'], r['tname'], str(score)]
	rows.append(row)

with open(args['out_tsv'], 'w') as out:
	for row in rows:
		out.write('\t'.join([str(_) for _ in row])+'\n')

