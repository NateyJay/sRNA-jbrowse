#!/usr/bin/env python3


import sys
import os
from collections import Counter
from pprint import pprint
from subprocess import call, Popen, PIPE


gff_input = sys.argv[1]
clean_input = gff_input.replace(".gff3", ".clean.gff3")
sorted_input = clean_input.replace(".gff3", ".sorted.gff3")
zipped_input = sorted_input.replace(".gff3", ".gff3.gz")


print(gff_input)
print(clean_input)
print(sorted_input)
print(zipped_input)


print("cleaning...")
c1 = f"python3 cleanup_gff.py -g {gff_input}"
call(c1.split())

print("sorting...")
c2 = f"gt gff3 -retainids -sortlines -tidy {clean_input}"
with open(sorted_input, 'w') as f:
	c2 = Popen(c2.split(), stdout=f)

print("zipping...")
c3 = f"bgzip {sorted_input}"
call(c3.split())

print("indexing...")
c4 = f"tabix -f -p gff {zipped_input}"
call(c4.split())


