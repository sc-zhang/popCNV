#!/usr/bin/env python
import os
import sys


def merge_CN(in_dir, out_cn):
	flist = os.listdir(in_dir)
	sp_list = []
	pos_list = []
	data_db = {}
	is_first = True
	
	print("Reading cn")
	for fn in flist:
		if fn.split('.')[-1].lower() != 'cn':
			continue
		print("\treading %s"%fn)
		sp = fn.split('.')[0]
		sp_list.append(sp)
		with open(os.path.join(in_dir, fn), 'r') as fin:
			for line in fin:
				data = line.strip().split()
				if is_first:
					pos_list.append(data[:3])
				key = "-".join(data[:3])
				if key not in data_db:
					data_db[key] = {}
				data_db[key][sp] = data[-1]
		is_first = False
	
	print("Writing data")
	with open(out_cn, 'w') as fout:
		fout.write("Chrom\tStart\tEnd\t%s\n"%('\t'.join(sp_list)))
		for pos in pos_list:
			fout.write("\t".join(pos))
			key = '-'.join(pos)
			for sp in sp_list:
				fout.write("\t"+data_db[key][sp])
			fout.write("\n")
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_dir> <out_cn>")
	else:
		in_dir, out_cn = sys.argv[1:]
		merge_CN(in_dir, out_cn)
