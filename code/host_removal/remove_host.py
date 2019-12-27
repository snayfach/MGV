
import csv, Bio.SeqIO

def id_host_region(gene2fam, num_genes, side):

	# init data
	block_size = 0
	gene_num = None
	num_pfams = 0
	
	# sort genes depending on start side
	if side == 'right':
		my_range = range(1,num_genes+1)[::-1]
	else:
		my_range = range(1,num_genes+1)
	
	# loop over genes
	for i in my_range:
		
		# unannotated gene: move to next
		if str(i) not in gene2fam:
			continue

		# viral gene: end of host region
		elif not 'PF' in gene2fam[str(i)]:
			return gene_num, num_pfams, block_size
	
		# non-viral gene: update host region boundary
		else:
			block_size += 1
			num_pfams += 1
			gene_num = i
	
	return gene_num, num_pfams, block_size


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('--in_base', required=True, metavar='PATH')
	parser.add_argument('--out_base', required=True, metavar='PATH')
	args = vars(parser.parse_args())
	
	# read gene coords
	gff = {}
	with open(in_gff) as f:
		next(f)
		contig_id = None
		for line in f:
			if line[0]=='#': continue
			contig_id, program, type, start, stop, gc, strand, phase, attr = line.rstrip('\n').split('\t')
			attr = dict([_.split('=') for _ in attr.rstrip(';').split(';')])
			gene_num = attr['ID'].split('_')[-1]
			r = {}
			r['contig_id'] = contig_id
			
			gff[contig_id+'_'+gene_num] = int(start), int(stop)

	# identify host regions
	info = {}
	for r in csv.DictReader(open(args['in_base']+'.tsv'), delimiter='\t'):
		
		info[r['contig_id']] = r
		
		# sequence boundaries
		fna_start, fna_stop = 1, int(r['length'])
		ffn_start, ffn_stop = 1, int(r['genes'])
		trim_left, trim_right = 'No', 'No'
		
		# fetch gene id of leftmost and rightmost non-viral gene
		# these determine host boundaries
		gene2fam = dict([_.split(':') for _ in r['gene2fam'].split(',')]) if r['gene2fam'] != '' else {}
		num_genes = int(r['genes'])

		# update left (5') boundary
		gene_id, num_pfams, num_genes = id_host_region(gene2fam, num_genes, 'left')
		if gene_id and num_pfams >= 2 and 100.0*num_pfams/num_genes >= 30.0:
			gene_start, gene_stop = coords[r['contig_id']+'_'+str(gene_id)]
			fna_start = gene_stop + 1
			ffn_start = gene_id + 1
			trim_left = 'Yes'
		
		# update right (3') boundary
		gene_id, num_pfams, num_genes = id_host_region(gene2fam, num_genes, 'right')
		if gene_id and num_pfams >= 2 and 100.0*num_pfams/num_genes >= 30.0:
			gene_start, gene_stop = coords[r['contig_id']+'_'+str(gene_id)]
			fna_stop = gene_start - 1
			ffn_stop = gene_id - 1
			trim_right = 'Yes'

		# update record
		info[r['contig_id']]['fna_start'] = fna_start
		info[r['contig_id']]['fna_stop'] = fna_stop
		info[r['contig_id']]['gene_start'] = ffn_start
		info[r['contig_id']]['gene_stop'] = ffn_stop
		info[r['contig_id']]['trim_left'] = trim_left
		info[r['contig_id']]['trim_right'] = trim_left


	# write trimmed fna
	with open(args['out_base']+'.fna', 'w') as f:
		for r in Bio.SeqIO.parse(args['in_base']+'.fna', 'fasta'):
			start = info[r.id]['fna_start'] - 1
			stop = info[r.id]['fna_stop']
			seq = str(r.seq[start:stop])
			f.write('>'+r.id+'\n'+seq+'\n')

	# write trimmed faa
	with open(args['out_base']+'.faa', 'w') as f:
		for r in Bio.SeqIO.parse(args['in_base']+'.faa', 'fasta'):
			contig_id, gene_num = r.id.rsplit('_', 1)
			gene_start = info[contig_id]['gene_start']
			gene_stop = info[contig_id]['gene_stop']
			if int(gene_num) >= gene_start and int(gene_num) <= gene_stop:
				f.write('>'+r.id+'\n'+str(seq)+'\n')

	# write trimmed ffn
	with open(args['out_base']+'.ffn', 'w') as f:
		for r in Bio.SeqIO.parse(args['in_base']+'.ffn', 'fasta'):
			contig_id, gene_num = r.id.rsplit('_', 1)
			gene_start = info[contig_id]['gene_start']
			gene_stop = info[contig_id]['gene_stop']
			if int(gene_num) >= gene_start and int(gene_num) <= gene_stop:
				f.write('>'+r.id+'\n'+str(seq)+'\n')

	# write trimmed gff
	# a bit tricky
	for r in coords:
		pass





