#!/usr/bin/env python
import sys
import numpy as np


def search_pos(pos_list, pos):
	s = 0
	e = len(pos_list)-1
	while s<=e:
		mid = (s+e)/2
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


def get_ctg_cn(in_cn, in_agp, out_cn):
	cn_db = {}
	print("Reading CN")
	with open(in_cn, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			sp = int(data[1])
			ep = int(data[2])
			cn = float(data[3])
			if chrn not in cn_db:
				cn_db[chrn] = []
			cn_db[chrn].append([sp, ep, cn])
	
	ctg_regions = {}
	print("Reading AGP")
	with open(in_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			if chrn[:3].lower() != 'chr' or data[4] == 'U':
				continue
			sp = int(data[1])
			ep = int(data[2])
			tig = data[5]
			ctg_regions[tig] = [chrn, sp, ep]
	
	print("Calculating CN")
	ctg_cn = {}
	for ctg in ctg_regions:
		chrn, sp, ep = ctg_regions[ctg]
		nsp = search_pos(cn_db[chrn], sp)
		nep = search_pos(cn_db[chrn], ep)
		if nsp == -1 or nep == -1:
			print("Error on %s\t%d\t%d"%(chrn, sp, ep))
			continue
		cn_list = []
		for i in range(nsp, nep+1):
			cn_list.append(cn_db[chrn][i][2])
		ctg_cn[ctg] = np.median(cn_list)
	
	print("Writing result")
	with open(out_cn, 'w') as fout:
		for ctg in sorted(ctg_cn):
			fout.write("%s\t%f\n"%(ctg, ctg_cn[ctg]))
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_cn> <in_agp> <out_cn>")
	else:
		in_cn, in_agp, out_cn = sys.argv[1:]
		get_ctg_cn(in_cn, in_agp, out_cn)
