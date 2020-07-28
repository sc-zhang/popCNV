#!/usr/bin/env python
import os
import argparse
import numpy as np
import scipy.stats


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument("-i", "--input", help="Input CN file", required=True)
	group.add_argument("-o", "--outdir", help="Output directory of results", required=True)
	group.add_argument("-g", "--group", help="Group id with populations", required=True)
	group.add_argument("-s", "--sample", help="Samples with groups", required=True)
	group.add_argument("-w", "--wild", help="Group id of wild population", required=True)
	return group.parse_args()


def calc_RFD(in_cn, grp_file, smp_file, wild_id, out_dir):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	print("Reading populations")
	pop_db = {}
	with open(grp_file, 'r') as fin:
		for line in fin:
			data = line.strip().split('\t')
			pop_db[data[0]] = data[1]
	
	print("Reading samples")
	sample_db = {}
	with open(smp_file, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			sample_db[data[0]] = data[1]
	
	
	print("Reading CN and calc RFD")
	RFD_db = {}
	sample_idx = {}
	with open(in_cn, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[0] == 'Gene':
				for i in range(1, len(data)):
					sample_idx[i] = data[i]
			else:
				gn = data[0]
				cn_db = {'pop': []}
				for i in range(1, len(data)):
					if data[i] == '-':
						continue
					val = int(data[i])
					sn = sample_idx[i]
					if sn not in sample_db:
						continue
					pid = sample_db[sn]
					if pid not in cn_db:
						cn_db[pid] = []
					cn_db[pid].append(val)
					cn_db['pop'].append(val)
				F_db = {}
				if len(cn_db['pop']) == 0:
					continue
				F_pop = max(np.bincount(cn_db['pop']))*1.0/len(cn_db['pop'])
				for pid in pop_db:
					if pid not in cn_db:
						continue
					if len(cn_db[pid]) == 0:
						continue
					F_db[pid] = max(np.bincount(cn_db[pid]))*1.0/len(cn_db[pid])
				
				if wild_id not in cn_db:
					continue
				for pid in pop_db:
					if pid == wild_id or pid not in cn_db:
						continue
					f, p = scipy.stats.f_oneway(cn_db[pid], cn_db[wild_id])

					RFD = (F_db[pid]-F_db[wild_id])*1.0/F_pop
					if pid not in RFD_db:
						RFD_db[pid] = []
					RFD_db[pid].append([RFD, p, gn])
	
	print("Getting gene list with max RFD")
	gene_RFD_db = {}
	for pid in RFD_db:
		cnt = 0
		with open(os.path.join(out_dir, pop_db[pid]+'.list'), 'w') as fout:
			for RFD, p, gn in sorted(RFD_db[pid], reverse=True):
				if gn not in gene_RFD_db:
					gene_RFD_db[gn] = {}
				gene_RFD_db[gn][pid] = [RFD, p]
				if cnt < len(RFD_db[pid])*0.05 and p<0.05:
					fout.write("%s\t%f\t%f\n"%(gn, RFD, p))
					cnt += 1
	print("Writing total RFDs")
	with open(os.path.join(out_dir, "Total.list"), 'w') as fout:
		fout.write("Gene")
		for pid in sorted(RFD_db):
			fout.write("\t%s"%(pop_db[pid]))
		fout.write("\n")
		
		for gene in sorted(gene_RFD_db):
			fout.write(gene)
			for pid in sorted(RFD_db):
				if pid in gene_RFD_db[gene]:
					fout.write("\t%f,%f"%(gene_RFD_db[gene][pid][0], gene_RFD_db[gene][pid][1]))					
				else:
					fout.write("\tnan,nan")
			fout.write("\n")
	print("Finished")


if __name__ == "__main__":
	opts = get_opts()
	in_cn = opts.input
	grp_file = opts.group
	smp_file = opts.sample
	wild_id = opts.wild
	out_dir = opts.outdir
	calc_RFD(in_cn, grp_file, smp_file, wild_id, out_dir)
