#!/usr/bin/env python3


import sys
import os
from collections import Counter
from pprint import pprint
# from math import floor
# from subprocess import Popen, PIPE
# from pathlib import Path

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gff_file', 
	required=True,
	help='gff3 file from shortstack')

args = parser.parse_args()
gff_file = args.gff_file



out_file = gff_file.replace(".gff3", ".clean.gff3")

len_d = {}

##sequence-region   TRIATcontig_18 1814 1415423

with open(gff_file, 'r') as f:
	lines = f.readlines()

header = lines.pop(0).strip()

for line in lines:

	line = line.strip().split('\t')

	chrom = line[0]
	start = int(line[3])
	stop  = int(line[4])

	try:
		len_d[chrom]
	except KeyError:
		len_d[chrom] = [start, stop]


	if len_d[chrom][0] > start:
		len_d[chrom][0] = start


	if len_d[chrom][1] < stop:
		len_d[chrom][1] = stop


for chrom, entry in len_d.items():
	start, stop = entry

	header += f"\n##sequence-region   {chrom} {start} {stop}"


with open(out_file, 'w') as outf:
	print(header, file=outf)

	for line in lines:

		line = line.strip()

		line = line.replace('DicerCall', 'dicerCall')
		line = line.replace('MIRNA', 'miRNA')


		print(line, file=outf)

# gt gff3 -retainids -sortlines -tidy $2 > $2

# bgzip $2
# tabix -f -p gff $1.gz 






