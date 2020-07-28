#!/usr/bin/env python
import sys


def rounding_cn(in_cn, out_cn):
	cnt = 0
	with open(in_cn, 'r') as fin:
		with open(out_cn, 'w') as fout:
			for line in fin:
				data = line.strip().split()
				if data[0] == 'Gene':
					fout.write(line)
				else:
					for i in range(1, len(data)):
						if data[i] == 'nan' or data[i] == 'inf':
							cnt += 1
							data[i] = '-'
						else:
							data[i] = str(int(round(float(data[i]), 0)))
					fout.write("\t".join(data)+"\n")
	print("Unvalid data count: %d"%cnt)


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_CN> <out_CN>")
	else:
		in_cn, out_cn = sys.argv[1:]
		rounding_cn(in_cn, out_cn)
