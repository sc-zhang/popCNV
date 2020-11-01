#!/usr/bin/env python
import sys
import os
import pysam
import multiprocessing


def calc_depth(genome_size, in_bam, sp):
	read_size = 0
	with pysam.AlignmentFile(in_bam, 'rb') as fin:
		for line in fin:
			ql = line.query_length
			if ql != -1:
				read_size += ql
	return read_size*1.0/genome_size, sp


def get_seq_depth(in_genome, in_bam_dir, out_depth, ts):
	print("Getting genome size")
	genome_size = 0
	with open(in_genome, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				continue
			else:
				genome_size += len(line.strip())
	
	print("Calculating seq depth")
	pool = multiprocessing.Pool(processes = ts)
	res = []
	for bam in os.listdir(in_bam_dir):
		if bam.split('.')[-1].lower() != 'bam':
			continue
		print("\tCalc %s"%bam)
		sp = bam.split('.')[0]
		r = pool.apply_async(calc_depth, args=(genome_size, os.path.join(in_bam_dir, bam), sp, ))
		res.append(r)
	pool.close()
	pool.join()
	
	seq_depth = []
	for r in res:
		depth, sp = r.get()
		seq_depth.append([sp, depth])
	
	print("Writing result")
	with open(out_depth, 'w') as fout:
		for sp, depth in sorted(seq_depth):
			fout.write("%s\t%f\n"%(sp, depth))

	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <in_genome> <in_bam_dir> <out_depth> <threads>")
	else:
		in_genome, in_bam_dir, out_depth, ts = sys.argv[1:]
		ts = int(ts)
		get_seq_depth(in_genome, in_bam_dir, out_depth, ts)
