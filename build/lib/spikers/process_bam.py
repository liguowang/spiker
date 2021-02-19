#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:00:53 2020

@author: m102324
"""
import pysam
import numpy
from scipy import stats


def bam_info (bamfile, layout, frac = 0.2, n=500000):
	'''
	Extract DNA fragment size, read length and chrom sizes information from
	BAM file. For PE, fragment size is estiamted from read pairs.
	For SE, cannot estimate fragment size from the BAM file, therefore,
	fragment size is set to None.

	Parameters
	----------
	bamfile : str
		Input BAM files.
	layout : str
		Must be "PE" (paired end) or "SE" (single end).
	n : TYPE, int
		Number of paired-end alignments sampled. The default is 500000.
	frac : float
		Fraction to cut off of both tails of the distribution.

	Returns
	-------
	int
		The median fragment size.

	'''
	samfile = pysam.AlignmentFile(bamfile,'rb')
	chrom_sizes = dict(zip(samfile.references, samfile.lengths))
	frag_counts = 0
	frag_sizes = []
	read_length = []
	if layout == 'PE':
		try:
			while (1):
				aligned_read = next(samfile)
				if aligned_read.is_qcfail:continue
				if aligned_read.is_duplicate:continue
				if aligned_read.is_secondary:continue
				if aligned_read.is_supplementary:continue
				if aligned_read.is_unmapped:continue
				if aligned_read.mate_is_unmapped:continue
				if aligned_read.is_read2:continue
				frag_counts +=1
				frag_sizes.append(abs(aligned_read.template_length))
				read_length.append(aligned_read.query_alignment_length)
				#print (aligned_read.query_name + '\t' + str(aligned_read.template_length) )
				if frag_counts > n:
					break
		except StopIteration:
			pass
		#the order of chroms must be consistent with those in bedGraph file.
		#chrom_sizes = sorted(list(zip(samfile.references, samfile.lengths)))
		return (stats.trim_mean(frag_sizes, frac), numpy.mean(read_length), chrom_sizes)
	elif layout == 'SE':
		try:
			while (1):
				aligned_read = next(samfile)
				if aligned_read.is_qcfail:continue
				if aligned_read.is_duplicate:continue
				if aligned_read.is_secondary:continue
				if aligned_read.is_supplementary:continue
				if aligned_read.is_unmapped:continue
				frag_counts +=1
				read_length.append(aligned_read.query_alignment_length)
				if frag_counts > n:
					break
		except StopIteration:
			pass
		return (None, numpy.mean(read_length), chrom_sizes)


def get_layout(bamfile, n=500000):
	'''
	Return sequencing layout

	Parameters
	----------
	bamfile : str
		Input BAM files.
	n : TYPE, int
		Number of paired-end alignments sampled. The default is 500000.

	Returns
	-------
	str
		'PE', 'SE', "Mix" or 'Unknown'

	'''
	samfile = pysam.AlignmentFile(bamfile,'rb')
	frag_counts = 0
	layout = {'pair':0, 'single':0}
	try:
		while (1):
			aligned_read = next(samfile)
			if aligned_read.is_qcfail:continue
			if aligned_read.is_duplicate:continue
			if aligned_read.is_secondary:continue
			if aligned_read.is_supplementary:continue
			if aligned_read.is_unmapped:continue
			frag_counts +=1
			if aligned_read.is_paired:
				layout['pair'] += 1
			else:
				layout['single'] += 1
			if frag_counts > n:
				break
	except StopIteration:
		print("Done")
	if layout['pair'] == frag_counts:
		return 'PE'
	elif layout['single'] == frag_counts:
		return 'SE'
	elif layout['pair'] > 0 and layout['single'] >0:
		return 'Mix'
	else:
		return 'Unknown'


def total_mapped(bamfile):
	"""
	Return total mapped reads.

	Parameters
	----------
	bamfile : str
		Input BAM files.

	Returns
	-------
	int
		The total mapped reads in BAM file.

	"""
	samfile = pysam.Samfile(bamfile,'rb')
	return samfile.mapped

if __name__=='__main__':
	import sys
	avg_d_size = bam_info(sys.argv[1], 100000)
	print (avg_d_size)