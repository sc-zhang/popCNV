#!/usr/bin/env python
import sys
import gzip
import argparse
import multiprocessing


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument('-g', '--genome', help='Input genome fasta file', required=True)
	group.add_argument('-b', help='window size or bed file contain windows', required=True)
	group.add_argument('-o', '--output', help='Output GC content file', required=True)
	group.add_argument('-t', '--threads', type=int, help='Threads', required=True, default=10)
	return group.parse_args()


def get_seq_GC(seq, bins, chrn):
	print("\tStating %s"%chrn)
	gc_list = []
	for sp, ep in bins:
		gc_cnt = 0
		for i in range(sp-1, ep):
			if seq[i].lower() == 'g' or seq[i].lower() == 'c':
				gc_cnt += 1
		gc_list.append(int(gc_cnt*100.0/(ep-sp+1)))
	#res_str = chrn+'\t'+'\t'.join(list(map(str, gc_list)))
	return gc_list, chrn


def GC_stat(genome_file, bed_or_size, out_file, ts):
	print("Reading genome")
	genome_db = {}
	with open(genome_file, 'r') as fin:
		seq = ''
		id = ''
		for line in fin:
			if line[0] == '>':
				if seq != '':
					genome_db[id] = seq
				id = line.strip()[1:]
				seq = ''
			else:
				seq += line.strip()
		genome_db[id] = seq
	
	print("Generating bins")
	bins = {}
	if bed_or_size.isdigit():
		bs = int(bed_or_size)
		for chrn in genome_db:
			bins[chrn] = []
			i = 1
			while i<=len(genome_db[chrn]):
				e = i+bs-1
				if e >= len(genome_db[chrn]):
					e = len(genome_db[chrn])
				bins[chrn].append([i, e])
				i += bs
	else:
		with open(bed_or_size, 'r') as fin:
			for line in fin:
				data = line.strip().split()
				chrn = data[0]
				sp = int(data[1])
				ep = int(data[2])
				if sp > ep:
					sp, ep = ep, sp
				bins[chrn].append([sp, ep])
			for chrn in bins:
				bins[chrn] = sorted(bins[chrn])
	
	print("Starting stat")
	pool = multiprocessing.Pool(processes = ts)
	res = []
	for chrn in genome_db:
		r = pool.apply_async(get_seq_GC, args=(genome_db[chrn], bins[chrn], chrn, ))
		res.append(r)
	pool.close()
	pool.join()
	
	gc_db = {}
	for r in res:
		gc_list, chrn = r.get()
		#data = res_str.split('\t')
		#chrn = data[0]
		#gc_list = list(map(int, data[1:]))
		gc_db[chrn] = gc_list
	
	print("Writing result")
	with open(out_file, 'w') as fout:
		for chrn in sorted(gc_db):
			for i in range(0, len(bins[chrn])):
				fout.write("%s\t%d\t%d\t%d\n"%(chrn, bins[chrn][i][0], bins[chrn][i][1], gc_db[chrn][i]))
	
	print("Finished")


if __name__ == "__main__":
	opts = get_opts()
	genome_file = opts.genome
	bed_or_size = opts.b
	out_file = opts.output
	ts = opts.threads
	GC_stat(genome_file, bed_or_size, out_file, ts)
