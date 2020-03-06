#!/usr/bin/env python
import os
import sys
import gzip
import multiprocessing



def sub_norm_rd(in_file, out_file, fn, sn, s_ratio):
	print("\tConverting\t%s"%fn)
	with gzip.open(in_file, 'rt') as fin:
		with open(out_file ,'w') as fout:
			for line in fin:
				data = line.strip().split()
				data[1] = str(int(data[1])+1)
				rd = float(data[-1])/s_ratio[sn]
				data[-1] = str(rd)
				fout.write("%s\n"%('\t'.join(data)))


def norm_rd(in_dir, sd, out_dir, ts):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	print("Getting sample ratio")
	s_ratio = {}
	with open(sd, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			s_ratio[data[0]] = float(data[1])
	
	print("Reading depth files and writing new depth files")
	pool = multiprocessing.Pool(processes = ts)
	res = []
	for fn in os.listdir(in_dir):
		if not fn.endswith(".regions.bed.gz"):
			continue
		sn = fn.split('.')[0]
		if sn not in s_ratio:
			continue
		in_file = os.path.join(in_dir, fn)
		out_file = os.path.join(out_dir, sn+".rd")
		r = pool.apply_async(sub_norm_rd, args=(in_file, out_file, fn, sn, s_ratio, ))
		res.append(r)

	pool.close()
	pool.join()
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <in_dir> <sample_ratio> <out_dir> <threads>")
	else:
		in_dir, sd, out_dir, ts = sys.argv[1:]
		ts = int(ts)
		norm_rd(in_dir, sd, out_dir, ts)
