#!/usr/bin/env python
import sys
import os
import numpy as np


def search_pos(pos_list, pos):
	s = 0
	e = len(pos_list)-1
	while s<=e:
		mid = int((s+e)/2)
		if pos_list[mid][0]>pos:
			e = mid-1
		elif pos_list[mid][0]<pos:
			s = mid+1
		else:
			return mid
	if pos_list[e][1] >= pos:
		return e
	else:
		return -1


def get_gene_cn(in_cn_dir, in_bed, out_cn):
	cn_db = {}
	sp_list = []
	print("Reading CN")
	for in_cn in os.listdir(in_cn_dir):
		if in_cn.split('.')[-1].lower() != 'cn':
			continue
		sa = in_cn.split('.')[0]
		cn_db[sa] = {}
		print("\tReading %s"%sa)
		sp_list.append(sa)
		with open(os.path.join(in_cn_dir, in_cn), 'r') as fin:
			for line in fin:
				data = line.strip().split()
				chrn = data[0]
				sp = int(data[1])
				ep = int(data[2])
				cn = float(data[3])
				if chrn not in cn_db[sa]:
					cn_db[sa][chrn] = []
				cn_db[sa][chrn].append([sp, ep, cn])

	sp_list = sorted(sp_list)
	gn_regions = {}
	print("Reading bed")
	with open(in_bed, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			sp = int(data[1])
			ep = int(data[2])
			gn = data[3]
			gn_regions[gn] = [chrn, sp, ep]
	
	print("Calculating CN")
	gn_cn = {}
	for gn in gn_regions:
		gn_cn[gn] = {}
		for sa in sp_list:
			chrn, sp, ep = gn_regions[gn]
			nsp = search_pos(cn_db[sa][chrn], sp)
			nep = search_pos(cn_db[sa][chrn], ep)
			if nsp == -1 or nep == -1:
				print("Error on %s\t%d\t%d"%(chrn, sp, ep))
				continue
			cn_list = []
			for i in range(nsp, nep+1):
				cn_list.append(cn_db[sa][chrn][i][2])
			try:
				gn_cn[gn][sa] = np.median(cn_list)
			except Exception:
				print(sa, cn_list)
	
	print("Writing result")
	with open(out_cn, 'w') as fout:
		fout.write("Gene\t%s\n"%('\t'.join(sp_list)))
		for gn in sorted(gn_cn):
			fout.write(gn)
			for sp in sp_list:
				if sp in gn_cn[gn]:
					fout.write("\t%f"%(gn_cn[gn][sp]))
				else:
					fout.write("\tnan")
			fout.write('\n')
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_cn_dir> <in_bed> <out_cn>")
	else:
		in_cn_dir, in_bed, out_cn = sys.argv[1:]
		get_gene_cn(in_cn_dir, in_bed, out_cn)
