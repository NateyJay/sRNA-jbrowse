#!/usr/bin/env python3


import sys
import os
from collections import Counter
# from pprint import pprint
# from math import floor
from subprocess import Popen, PIPE
from pathlib import Path

# from multiprocessing import Lock, Process, Queue, current_process, Pool
# import time
import shutil
# from datetime import datetime
from os.path import isfile, isdir

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam_file', 
	nargs="?",
	required=True,
	help='bamfile of aligned small RNAs (tested with shortstack)')

parser.add_argument('-r', '--readgroups', 
	nargs='+', 
	help='list of read_groups to be read. Specify `all` to report all. Check RGs with `samtools view -H bamfile.bam`', 
	default = [False])

parser.add_argument('-c', '--chromosomes', 
	nargs='+', 
	help='list of chromosomes to be read. Defaults to report all. Check chromosomes with `samtools view -H bamfile.bam`.', 
	default = [False])

parser.add_argument('-s', '--sizes', 
	nargs='+', 
	help='list of sRNA sizes to be analyzed separately', 
	default = [20,21,22,23,24])

parser.add_argument('-o', '--output_directory', 
	nargs="?",
	required=True,
	help='folder to save resulting wigs/bigwigs')

parser.add_argument("--delete_wigs",
	action='store_true',
	default=False,
	help='include to delete wig files from the result')

parser.add_argument("--quiet",
	action='store_true',
	default=False,
	help='prints output in a much more discrete format')



args = parser.parse_args()
input_bam  = args.bam_file
out_dir    = args.output_directory
rgs        = args.readgroups
sizes      = args.sizes
del_wigs   = args.delete_wigs
input_chroms = args.chromosomes
quiet      = args.quiet

sizes = [int(s) for s in sizes]
sizes = [s for s in range(min(sizes), max(sizes)+1)]



print("")
print(f"bigwigify.py")
print()
print(f"       input_bam: {input_bam}")
print(f"output_directory: {out_dir}")
print(f"      readgroups: {rgs}")
print(f"     chromosomes: {input_chroms}")
print(f"           sizes: {sizes}")
print(f"     delete wigs: {del_wigs}", flush=True)


# input_bam = "../01out-alignment/merged_alignments.bam"
# out_dir = 'bigwig'
# sizes = [20,21,22,23,24]

# rgs =  ['B-A1.t',
#  'B-A2.t',
#  'B-A3.t']


out_dir = out_dir.rstrip("/")
assert not isdir(out_dir), "Output dir exists! (will not overwrite)"




def get_chromosomes(file):
	chromosomes = []
	rgs = []
	call = f"samtools view -@ 4 -H {file}"
	# print(call)

	p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')
	out, err = p.communicate()

	for o in out.strip().split("\n"):
		o = o.strip().split('\t')
		if o[0] == "@SQ":
			name = o[1].split(":")[-1]
			length = int(o[2].split(":")[-1])
			chromosomes.append((name,length))
		if o[0] == "@RG":
			rgs.append(o[1].split(":")[-1])

	return(chromosomes, rgs)



chromosomes, bam_rgs = get_chromosomes(input_bam)

if input_chroms[0]:
	chromosomes = [(c,l) for c,l in chromosomes if c in input_chroms]

for c in input_chroms:
	if c == False:
		input_chroms = [ch for ch,l in chromosomes]
		break
	if c not in [ch for ch,l in chromosomes]:
		sys.exit(f"\nERROR: input chromosome not listed in bamfile\n{c}\n{chromosomes}")

# pprint(chromosomes)
# print(input_chroms)

if rgs[0] == 'all':
	rgs = bam_rgs
elif rgs[0] == False:
	print("\nERROR: readgroups not specified...")
	print("must come from this list:")
	print(bam_rgs)
	sys.exit()



for rg in rgs:
	assert rg in bam_rgs, f"{rg} not found in bamfile...\n{bam_rgs}"


Path(out_dir).mkdir(parents=True, exist_ok=True)

with open("chrom.sizes.txt", 'w') as outf:
	for chrom, size in chromosomes:
		print(chrom, size, sep='\t', file=outf)





class wiggleClass():
	def __init__(self, strand, size):
		self.buffer = []
		self.start_pos = 1
		self.last_dep = 0
		self.file = f"{out_dir}/{strand}{size}.wig"

	def extend_buffer(self, length):

		for l in range(length):
			try:
				self.buffer[l] += 1
			except IndexError:
				self.buffer.append(1)


	def get_report(self, pos):

		buffer    = self.buffer
		last_dep  = self.last_dep
		start_pos = self.start_pos
		file      = self.file

		try:
			dep = buffer.pop(0)
		except IndexError:
			dep = 0


		if dep != last_dep:

			span = pos - start_pos

			if span > 0:

				self.last_dep  = dep
				self.start_pos = pos

				# print(file, span, start_pos, last_dep)

				return(file, span, start_pos, last_dep)

		else:
			# print(False)
			return(False)



def process_line(line):

	flag = line[1]

	if flag == "16":
		strand = '-'
	elif flag == '0':
		strand = "+"
	else:
		strand = False

	length = len(line[9])

	if length < min(sizes):
		size = "short"
	elif length > max(sizes):
		size = "long"
	else:
		size = length


	sam_pos = int(line[3])
	sam_chrom = line[2]

	return(flag, size, strand, length, sam_pos, sam_chrom)


def samtools_view(bam, rgs, chroms):

	for chrom in chroms:
		rg_opt = " -r ".join(rgs)
		call = f"samtools view -@ 4 -F 4 -r {rg_opt} {bam} {chrom}"
		# print(call)
		p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')


		for i,line in enumerate(p.stdout):

			# if i > 1000000:
			# 	return

			line = line.strip().split()

			flag, size, strand, length, sam_pos, sam_chrom = process_line(line)

			if strand:
				yield(flag, size, strand, length, sam_pos, sam_chrom)

		p.wait()



samtools_iter = samtools_view(input_bam, rgs, [c for c,l in chromosomes])

strands = ["+","-"]

file_names = []


depth_c = Counter()

wig_d = {}

keys = []
for size in sizes + ['short', 'long']:
	for strand in strands:
		key = (size, strand)
		keys.append(key)
		with open(f"{out_dir}/{size}{strand}.wig", 'w') as outf:
			outf.write("")

def wigsZeroed():
	for key in keys:
		if wig_d[key].last_dep > 0:
			# print(wig_d[key].buffer)
			# input()
			return(False)
	return(True)


depth = 0
flag = False
sam_empty = False

for chrom, chrom_length in chromosomes:

	# pos_iter = [*range(1, chrom_length+1)]

	if not quiet:
		print()
		print(chrom)
		print("  ", end='')

	for key in keys:
		size, strand = key

		# pos_d[key] = Counter()
		wig_d[key] = wiggleClass(size, strand)


	if not flag:
		flag, size, strand, length, sam_pos, sam_chrom = next(samtools_iter)
		depth += 1


	pos = 1

	loop = True

	while loop:

		if not quiet:
			if depth % 10000 == 0:
				print(".", end = '', flush=True)

		if sam_pos == pos and sam_chrom == chrom and not sam_empty:


			wig_d[(size,strand)].extend_buffer(length)

			try:
				flag, size, strand, length, sam_pos, sam_chrom = next(samtools_iter)
				depth += 1
			except StopIteration:
				sam_empty =True

		elif sam_pos > pos or chrom != sam_chrom or sam_empty:
			# print('herere')
			for key in keys:


				report = wig_d[key].get_report(pos)

				if report:
					file, span, start_pos, last_dep = report
					# print()
					# print(key)
					# print(file)
					# print(f"variableStep chrom={chrom} span={span}")
					# print(f"{start_pos} {last_dep}")
					# input()

					with open(file, 'a') as outf:

						print(f"variableStep chrom={chrom} span={span}", file=outf)
						print(f"{start_pos} {last_dep}", file=outf)



			zeroed = wigsZeroed()
			if zeroed:
				# print('buffer empty!')
				pos = sam_pos
			else:
				pos += 1

			# pos += 1

			# print('here')


			if chrom != sam_chrom or sam_empty:
				loop = False


				for key in keys:
					start_pos = wig_d[key].start_pos

					span = chrom_length - start_pos

					if span > 0:

						with open(wig_d[key].file, 'a') as outf:

							print(f"variableStep chrom={chrom} span={span}", file=outf)
							print(f"{start_pos} 0", file=outf)
		# print(chrom, pos, sam_empty, loop)





print()
print(f"overall depth: {depth:,}")

print()
print('calculating RPM and converting to bigwigs...')

for key in keys:
	size, strand = key

	wig        = f"{out_dir}/{size}{strand}.wig"
	scaled_wig = f"{out_dir}/{size}{strand}.scaled.wig"
	bigwig     = wig.replace(".wig", ".bigwig")

	with open(wig, 'r') as f:
		with open(scaled_wig, 'w') as outf:
			for line in f:
				if "variableStep" in line:
					print(line.strip(), file=outf)
				else:
					line = line.strip().split()
					line = [int(l) for l in line]
					pos, dep = line

					dep = round(float(dep) / depth * 1000000,3)

					if strand == '-':
						dep = dep * -1

					print(pos, dep, file=outf)


	os.remove(wig)
	os.rename(scaled_wig, wig)


	print(f"  {wig} -> {bigwig}", flush=True)


	call = f"wigToBigWig {wig} ./chrom.sizes.txt {bigwig}"

	p = Popen(call.split(), stdout=PIPE, stderr=PIPE, encoding='utf-8')

	out, err= p.communicate()

	if out.strip() + err.strip() != "":

		print(out)
		print(err)

	if del_wigs:
		os.remove(wig)













