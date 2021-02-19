#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:39:53 2021

@author: m102324
"""
import os
import deeptools.countReadsPerBin
import pysam

def frip_score(peak_file, bam_file, nthreads = 8):
	if (os.path.exists(peak_file) and os.path.getsize(peak_file) > 0):
		peak_files = [peak_file]
		bam_files = [bam_file]
		cr = deeptools.countReadsPerBin.CountReadsPerBin(bam_files, bedFile=peak_files, numberOfProcessors=nthreads)
		reads_at_peaks = cr.run()
		reads_in_peaks = reads_at_peaks.sum(axis=0)[0]

		tmp = pysam.AlignmentFile(bam_file)
		total_reads = tmp.mapped
		frip = float(reads_in_peaks)/total_reads
		return frip
	else:
		return 0