
import Bio.SeqIO, time, gzip
from operator import itemgetter
import numpy as np, argparse

def parse_blast(handle):
	for line in handle:
		r = line.split()
		yield {
			'qname':r[0],
			'tname':r[1],
			'pid':float(r[2]),
			'len':float(r[3]),
			'qcoords':sorted([int(r[6]), int(r[7])]),
			'tcoords':sorted([int(r[8]), int(r[9])]),
			'qlen':float(r[-2]),
			'tlen':float(r[-1]),
			'evalue':float(r[-4])
			}

def yield_alignment_blocks(handle):
	# init block with 1st record
	key, alns = None, None
	for aln in parse_blast(handle):
		if aln['qname'] == aln['tname']:
			continue
		key = (aln['qname'], aln['tname'])
		alns = [aln]
		break
	# loop over remaining records
	for aln in parse_blast(handle):
		# skip self hits
		if aln['qname'] == aln['tname']:
			continue
		# extend block
		elif (aln['qname'], aln['tname']) == key:
			alns.append(aln)
		# yield block and start new one
		else:
			yield alns
			key = (aln['qname'], aln['tname'])
			alns = [aln]
	yield alns

def prune_alns(alns, min_length=0, min_evalue=1e-3):
	# remove short aligns
	alns = [aln for aln in alns if aln['len'] >= min_length and aln['evalue'] <= min_evalue]
	return alns

def compute_ani(alns):
	return round(sum(a['len'] * a['pid'] for a in alns)/sum(a['len'] for a in alns),2)

def compute_cov(alns):

	# merge qcoords
	coords = sorted([a['qcoords'] for a in alns])
	nr_coords = [coords[0]]
	for start, stop in coords[1:]:
		# overlapping, update start coord
		if start <= (nr_coords[-1][1] + 1):
			nr_coords[-1][1] = max(nr_coords[-1][1], stop)
		# non-overlapping, append to list
		else:
			nr_coords.append([start, stop])
	# compute qry_cov
	alen = sum([stop - start + 1 for start, stop in nr_coords])
	qcov = round(100.0*alen/alns[0]['qlen'],2)

	# merge tcoords
	coords = sorted([a['tcoords'] for a in alns])
	nr_coords = [coords[0]]
	for start, stop in coords[1:]:
		# overlapping, update start coord
		if start <= (nr_coords[-1][1] + 1):
			nr_coords[-1][1] = max(nr_coords[-1][1], stop)
		# non-overlapping, append to list
		else:
			nr_coords.append([start, stop])
	# compute qry_cov
	alen = sum([stop - start + 1 for start, stop in nr_coords])
	tcov = round(100.0*alen/ alns[0]['tlen'],2)

	return qcov, tcov


def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', dest='input', type=str, required=True, metavar='PATH',
		help="path to blastn input file (format: 'std 6 qlen slen')")
	parser.add_argument('-o', dest='output', type=str, required=True, metavar='PATH',
		help="path to ani file")
	parser.add_argument('-l', dest='length', type=int, metavar='INT',
		help="minimum alignment length to keep")
	return vars(parser.parse_args())
		
if __name__ == "__main__":
	args = parse_arguments()
	out = gzip.open(args['output'], 'w') if args['output'].split('.')[-1]=='gz' else open(args['output'], 'w')
	fields = ['qname', 'tname', 'num_alns', 'pid', 'qcov', 'tcov']
	out.write('\t'.join(fields)+'\n')
	input = gzip.open(args['input']) if args['input'].split('.')[-1]=='gz' else open(args['input'])
	for alns in yield_alignment_blocks(input):
		alns = prune_alns(alns)
		if len(alns) == 0: continue
		qname, tname = alns[0]['qname'], alns[0]['tname']
		ani = compute_ani(alns)
		qcov, tcov = compute_cov(alns)
		row = [qname, tname, len(alns), ani, qcov, tcov]
		out.write('\t'.join([str(_) for _ in row])+'\n')

