#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:34:13 2021

@author: m102324
"""
import pyBigWig
import logging
import os

def bdg2bw(in_bdg, out_bw, chromSizes):
	"""
	Convert bedGraph format file into bigwig.

	Parameters
	----------
	in_bdg : str
		bedGraph format file.
	out_bw : str
		Output bigwig file.
	chromSizes : dict
		Dictionary of chromosome IDs and sizes.

	Returns
	-------
	None.

	"""
	if os.path.exists(out_bw) and os.path.getsize(out_bw) > 0:
		logging.warning("\t\"%s\" exists and non-empy, skip." % out_bw)
	else:
		# get the order of chrom ids in bdg file.
		chrom_list = []
		BDG = open(in_bdg,'r')
		for l in BDG:
			l = l.strip()
			if l.startswith('track'):
				continue
			chrom_id = l.split()[0]
			if chrom_id not in chrom_list:
				chrom_list.append(chrom_id)
			else:
				continue
		BDG.close()

		# create chrom_sizes (used as bigwig header). eg: [("chr1", 1000000), ("chr2", 1500000)]
		chrom_sizes = []
		for chrom in chrom_list:
			if chrom in chromSizes:
				chrom_sizes.append((chrom, chromSizes[chrom]))
			else:
				continue

		BWOUT = pyBigWig.open(out_bw, "w")
		# add header to bigwig file
		BWOUT.addHeader(chrom_sizes)
		#logging.info("Sort %s ..." % in_bdg)
		#subprocess.call("tail -n +2 %s | sort -k1,1 -k2,2n > %s" % (in_bdg, in_bdg + '.sorted'),shell=True)

		for l in open(in_bdg, 'r'):
			l = l.strip()
			if l.startswith('track'):
				continue
			else:
				chrom,start,end,score = l.split()
				BWOUT.addEntries([chrom], [int(start)], [int(end)], [float(score)])
		BWOUT.close()