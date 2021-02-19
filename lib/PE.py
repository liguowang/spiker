#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:18:56 2020

@author: m102324
"""
import os
import logging

def pe2se (in_PE, out_SE, rl):
	'''
	Convert paired-end BED into single-end BED.

	Parameters
	----------
	in_PE : str
		Input PE BED file. Must have at least 3 columns including "chrom", "start", "end".
	out_SE : str
		Output SE BED file.
	fl : int
		Read length.
	'''

	if os.path.exists(out_SE) and os.path.getsize(out_SE) > 0:
		logging.warning("\tThe SE bed file exists and non-empty: %s. Skip." % out_SE)
		pass
	else:
		OUT = open(out_SE,'w')
		line_count = 0
		with open(in_PE,'r') as F:
			for l in F:
				l = l.strip()
				f = l.split()
				if len(f) != 3:
					continue
				try:
					chrom = f[0]
					start = int(f[1])
					end = int(f[2])
				except:
					continue
				line_count += 1
				print ('\t'.join([str(i) for i in (chrom, start, int(start + rl), '.','.','+')]), file=OUT)
				if end - rl > 0:
					print ('\t'.join([str(i) for i in (chrom, int(end - rl), end, '.','.','-')]), file=OUT)
				else:
					print ('\t'.join([str(i) for i in (chrom, 0, end, '.','.','-')]), file=OUT)

if __name__=='__main__':
	import sys
	pe2se(in_PE = sys.argv[1],out_SE=sys.argv[2], rl= 51)
