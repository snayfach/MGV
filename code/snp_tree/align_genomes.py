#!/usr/bin/env python

import subprocess as sp, shutil, os, time, argparse, sys, random

def fetch_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--genomes', type=str, required=True,
		help="""Directory of genomes OR list of file paths to genomes""")
	parser.add_argument('--ref', type=str,
		help="""Path to reference genome (default=pick random genome)""")
	parser.add_argument('--out', type=str, required=True,
		help="""Path to output directory""")
	parser.add_argument('--max', type=int, default=None,
		help="""Maximum # of genomes to process; useful for quick tests""")
	return vars(parser.parse_args())

def check_exe():
	programs = ['nucmer', 'show-snps', 'show-coords', 'show-diff']
	for prog in programs:
		if not any(os.access(os.path.join(path, prog), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
			sys.exit("Error: program '%s' not found on $PATH" % prog)

if __name__ == "__main__":

	start = time.time()
	
	args = fetch_args()
		
	check_exe()
	
	# setup output directory
	if not os.path.exists(args['out']) and not args['overwrite']:
		os.makedirs(args['out'])

	# setup input paths
	print("Reading input...")

	# fetch paths
	paths = []
	if os.path.isdir(args['genomes']):
		for _ in os.listdir(args['genomes']):
			path = os.path.join(args['genomes'], _)
			if path not in paths:
				paths.append(path)
	else:
		for _ in open(args['genomes']):
			path = _.rstrip()
			if path not in paths:
				paths.append(path)

	# add ref; make it 1st in path list
	if args['ref']:
		if args['ref'] in paths:
			paths.remove(args['ref'])
	else:
		args['ref'] = random.sample(paths, 1)[0]
		paths.remove(args['ref'])
	paths = [args['ref']]+paths

	# subset paths
	if args['max']:
		paths = paths[0:args['max']]

	# check paths
	for path in paths:
		if not os.path.exists(path):
			sys.exit("Error: file does not exist: %s" % path)

	# check names
	if len(set([os.path.basename(_) for _ in paths])) != len(paths):
		sys.exit("Error: genome names are not unique: %s" % path)

	# copy ref
	shutil.copy(args['ref'], os.path.join(args['out'], 'reference.fna'))

	print("   total genomes: %s" % len(paths))
	print("   ref genome: %s" % os.path.basename(args['ref']))
	print("   output: %s" % args['out'])

	# main
	print("\nAligning genomes...")
	for path in paths :
		
		genome_id = path.split('/')[-1].replace('.fna', '')
		print("   %s" % genome_id)
		
		out_dir = '%s/aln/%s' % (args['out'], genome_id)
		if not os.path.exists(out_dir): os.makedirs(out_dir)
		
		log = open(out_dir+'/log','w')

		command = "nucmer %s " % args['ref']
		command += "%s " % path
		command += "--prefix %s/%s " % (out_dir, genome_id)
		process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
		out, err = process.communicate()
		log.write(str(out)+'\n'+str(err))
		
		for utility in ['coords', 'snps', 'diff']:
			command = "show-%s " % utility
			command += "%s/%s.delta " % (out_dir, genome_id)
			command += "> %s/%s" % (out_dir, utility)
			process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
			out, err = process.communicate()
			log.write(str(out)+'\n'+str(err))

	print("\nDone!")


