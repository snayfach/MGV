#!/usr/bin/env python

import os, sys

class ArrayCluster:
	def __init__(self, array):
		self.arrays = [array]
		self.contig_id = array.info['contig_id']
		self.start_pos = array.info['start_pos']
		self.end_pos = array.info['end_pos']
		
	def append(self, array):
		self.arrays.append(array)
		self.end_pos = array.info['end_pos']
		
	def sort_arrays(self):
		""" Sort arrays in cluster by array length (asc) and repeat length (asc) """
		coords = [(a.info['array_length'], a.info['repeat_length']) for a in self.arrays]
		indexes = sorted(range(len(coords)), key=lambda k: coords[k])
		self.arrays = [self.arrays[i] for i in indexes]
	
	def pick_arrays(self):
		""" Pick non-overlapping set of arrays from cluster """
		if len(self.arrays) == 1:
			chosen_array = self.arrays[0]
			yield chosen_array
		else:
			# sort arrays by array length and repeat length (asc)
			self.sort_arrays()
			while len(self.arrays) > 0:
				# choose last array
				chosen_array = self.arrays.pop() 
				# delete intersecting arrays
				for array_index, array in enumerate(self.arrays):
					if (array.info['start_pos'] >= chosen_array.info['start_pos'] and
							array.info['start_pos'] <= chosen_array.info['end_pos']):
						del self.arrays[array_index]
					elif (array.info['end_pos'] >= chosen_array.info['start_pos'] and
							array.info['end_pos'] <= chosen_array.info['end_pos']):	
						del self.arrays[array_index]	
				yield chosen_array
				
	def print_coords(self):
		print("\n\nCoordinates")
		print("#################################")
		for array in self.arrays:
			print(taxon_id, array.info['contig_id'], array.info['start_pos'], array.info['end_pos'], array.info['tool'])

	def print_spacers(self):
		print("Repeat/Spacers")
		print("#################################")
		for array_index, array in enumerate(self.arrays):
			print(taxon_id, array.info['contig_id'], array.info['start_pos'], array.info['end_pos'], array.info['tool'])
			for spacer in array.spacers:
				print(spacer['start_pos'], spacer['repeat_seq'], spacer['spacer_seq'])

class Array:
	def __init__(self, info, tool):
		self.info = info
		self.spacers = []
		self.info['start_pos'] = int(self.info['start_pos'])
		self.info['end_pos'] = int(self.info['end_pos'])
		self.info['array_length'] = self.info['end_pos'] - self.info['start_pos'] + 1
		self.info['repeat_length'] = len(self.info['consensus_repeat'])
		self.info['tool'] = tool

def read_arrays(basename, tool):
	import csv
	arrays = {}
	arrays_path = basename+'.arrays'
	spacers_path = basename+'.spacers'
	if (not os.path.exists(arrays_path) or 
			not os.path.exists(spacers_path)
		):
		sys.stderr.write("Warning: arrays missing: %s\n" % basename)
		return arrays
	with open(basename+'.arrays') as f:
		for d in csv.DictReader(f, delimiter='\t'):
			arrays[d['contig_id'], d['array_num']] = Array(d, tool)				
	with open(basename+'.spacers') as f:
		for d in csv.DictReader(f, delimiter='\t'):
			arrays[d['contig_id'], d['array_num']].spacers.append(d)
	return arrays.values()
		
def write_arrays(arrays, basename):
	if not os.path.isdir(os.path.dirname(basename)):
		os.makedirs(os.path.dirname(basename))

	arrays_file = open(basename+'.arrays', 'w')
	arrays_fields = ['contig_id', 'array_num', 'start_pos', 
					 'end_pos', 'truncated', 'count_spacers', 
			         'percent_id', 'repeat_length', 'consensus_repeat',
					 'tool']
	arrays_file.write('\t'.join(arrays_fields)+'\n')
	for array_num, array in enumerate(arrays):
		array.info['array_num'] = array_num + 1
		values = [str(array.info[f]) for f in arrays_fields]
		arrays_file.write('\t'.join(values)+'\n')
			
	spacers_file = open(basename+'.spacers', 'w')
	spacer_fields = ['contig_id', 'array_num', 'repeat_num', 
	                 'start_pos', 'repeat_seq', 'spacer_seq']
	spacers_file.write('\t'.join(spacer_fields)+'\n')
	for array_num, array in enumerate(arrays):
		for spacer in array.spacers:
			spacer['array_num'] = array_num + 1
			values = [str(spacer[f]) for f in spacer_fields]
			spacers_file.write('\t'.join(values)+'\n')

def sort_arrays(arrays):
	""" Sort arrays by contig_id (asc), start_pos (asc), end_pos (desc)"""
	coords = [(a.info['contig_id'], a.info['start_pos'], -a.info['end_pos']) for a in arrays]
	indexes = sorted(range(len(coords)), key=lambda k: coords[k])
	return [arrays[i] for i in indexes]
	
def fetch_taxon_ids():
	dir = '/global/projectb/scratch/snayfach/projects/pbin/0_input_data/0_select_taxa'
	ids = [(_.rstrip(), 'meta') for _ in open(dir+'/taxids_meta.txt')]
	ids += [(_.rstrip(), 'isolate') for _ in open(dir+'/taxids_isolate.txt')]
	return ids

class SummaryStats:
	def __init__(self):
		
		self.stats = {}
		for level in ['total', 'cluster', 'single', 'nr']:
			self.stats[level] = {}
			self.stats[level]['count'] = 0 
			for object in ['arrays', 'spacers']:
				self.stats[level][object] = {}
				for tool in ['crt', 'pilercr']:
					self.stats[level][object][tool] = 0
		
	def update_total(self, arrays):
		for array in arrays:
			self.stats['total']['count'] += 1
			self.stats['total']['arrays'][array.info['tool']] += 1
			self.stats['total']['spacers'][array.info['tool']] += len(array.spacers)-1
			
	def update_nr(self, arrays):
		for array in arrays:
			self.stats['nr']['count'] += 1
			self.stats['nr']['arrays'][array.info['tool']] += 1
			self.stats['nr']['spacers'][array.info['tool']] += len(array.spacers)-1
			
	def update_clusters(self, clusters):
		for cluster in clusters:
			type = 'single' if len(cluster.arrays) == 1 else 'cluster'
			self.stats[type]['count'] += 1
			for array in cluster.arrays:
				self.stats[type]['arrays'][array.info['tool']] += 1
				self.stats[type]['spacers'][array.info['tool']] += len(array.spacers)-1
			
	def print_report(self):

		total_arrays = sum(self.stats['total']['arrays'].values())
		total_spacers = sum(self.stats['total']['spacers'].values())
		print("Total arrays/spacers: %s/%s" % (total_arrays, total_spacers))

		crt_arrays = self.stats['total']['arrays']['crt']
		crt_spacers = self.stats['total']['spacers']['crt']
		print("  CRT: %s/%s" % (crt_arrays, crt_spacers))
		
		pilercr_arrays = self.stats['total']['arrays']['pilercr']
		pilercr_spacers = self.stats['total']['spacers']['pilercr']
		print("  PILERCR: %s/%s" % (pilercr_arrays, pilercr_spacers))

		print("\nSingletons/clusters: %s/%s" % (self.stats['single']['count'], self.stats['cluster']['count']))

		total_arrays = sum(self.stats['nr']['arrays'].values())
		total_spacers = sum(self.stats['nr']['spacers'].values())
		print("\nNR arrays/spacers: %s/%s" % (total_arrays, total_spacers))

		crt_arrays = self.stats['nr']['arrays']['crt']
		crt_spacers = self.stats['nr']['spacers']['crt']
		print("  CRT: %s/%s" % (crt_arrays, crt_spacers))
		
		pilercr_arrays = self.stats['nr']['arrays']['pilercr']
		pilercr_spacers = self.stats['nr']['spacers']['pilercr']
		print("  PILERCR: %s/%s" % (pilercr_arrays, pilercr_spacers))		
		

if __name__ == "__main__":
	
	stats = SummaryStats()
	
	crt_base = sys.argv[1]
	plr_base = sys.argv[2]
	out_base = sys.argv[3]
	
	# read in arrays from file
	arrays = []
	arrays += read_arrays(crt_base, 'crt')
	arrays += read_arrays(plr_base, 'pilercr')
	stats.update_total(arrays)
	
	# build array clusters
	clusters = []
	for index, array in enumerate(sort_arrays(arrays)):
		# start new array cluster
		if (index == 0 or # first array in list
				array.info['contig_id'] != clusters[-1].contig_id or # different contig
				array.info['start_pos'] > clusters[-1].end_pos # downstream
			):
			clusters.append(ArrayCluster(array))
		# add array to existing cluster and update coords
		else:
			clusters[-1].append(array)
	stats.update_clusters(clusters)
	
	# pick non-overlapping set of arrays per cluster
	nr_arrays = []
	for cluster in clusters:
		for array in cluster.pick_arrays():
			nr_arrays.append(array)
	stats.update_nr(nr_arrays)

	# one last sort
	sorted_arrays = sort_arrays(nr_arrays)

	# write non-overlapping arrays to file
	write_arrays(sorted_arrays, out_base)
		
	stats.print_report()

				

	
			

		 
