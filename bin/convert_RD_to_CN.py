#!/usr/bin/env python
import sys
import numpy as np


def convert_RD_to_CN(in_gc, in_rd, out_cn):
	info_db = {}
	print("Reading GC")
	with open(in_gc, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			key = '-'.join(data[:3])
			if key not in info_db:
				info_db[key] = {}
			info_db[key]['gc'] = int(data[-1])
	
	print("Reading RD")
	with open(in_rd, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			key = '-'.join(data[:3])
			info_db[key]['rd'] = float(data[-1])
	
	print("Calculating CN")
	gc_db = {}
	for key in info_db:
		GC = info_db[key]['gc']
		if GC not in gc_db:
			gc_db[GC] = []
		gc_db[GC].append(info_db[key]['rd'])
	
	cn_db = {}
	for GC in gc_db:
		cn = np.median(gc_db[GC])
		cn_db[GC] = cn

	print("Converting data")
	data_db = {}
	for key in info_db:
		chrn, sp, ep = key.split('-')
		if chrn not in data_db:
			data_db[chrn] = []
		sp = int(sp)
		ep = int(ep)
		cn = info_db[key]['rd']*1.0/cn_db[info_db[key]['gc']]
		data_db[chrn].append([sp, ep, cn])
	
	print("Writing result")
	with open(out_cn, 'w') as fout:
		for chrn in sorted(data_db):
			for sp, ep, cn in sorted(data_db[chrn]):
				fout.write("%s\t%d\t%d\t%f\n"%(chrn, sp, ep, cn))
	
	print("Finished")

if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_gc> <in_rd> <out_cn>")
	else:
		in_gc, in_rd, out_cn = sys.argv[1:]
		convert_RD_to_CN(in_gc, in_rd, out_cn)
