#!/usr/bin/env python
import os
import sys


def norm_rd(in_dir, sd, out_dir):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	print("Getting sample ratio")
	s_ratio = {}
	with open(sd, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			s_ratio[data[0]] = float(data[1])
	
	print("Reading depth files and writing new depth files")
	for fn in os.listdir(in_dir):
		sn = fn.split('.')[0]
		if sn not in s_ratio:
			continue
		with open(os.path.join(in_dir, fn), 'r') as fin:
			with open(os.path.join(out_dir, fn) ,'w') as fout:
				print("\tConverting\t%s"%fn)
				for line in fin:
					data = line.strip().split()
					rd = float(data[-1])/s_ratio[sn]
					data[-1] = str(rd)
					fout.write("%s\n"%('\t'.join(data)))	
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_dir> <sample_ratio> <out_dir>")
	else:
		in_dir, sd, out_dir = sys.argv[1:]
		norm_rd(in_dir, sd, out_dir)
