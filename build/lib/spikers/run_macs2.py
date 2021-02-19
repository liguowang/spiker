#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 15:59:28 2020
@author: m102324
"""
import sys,os
import logging
import shutil
import subprocess
from numpy import log10
import random
import string

def find_macs2():
	'''
	Find the path of macs2 command.
	'''
	cmd = shutil.which("macs2")
	if cmd is None:
		logging.error("\tCannot find the \"macs2\" command!")
		sys.exit()
	logging.debug("\tFound the macs2 command:\"%s\"" %  cmd)
	return cmd

def deduplicate(bamfile, outfile, informat,  gsize='hs', keepdup = 1, verbose=0):
	'''
	Remove the duplicate reads from BAM file.

	parameters
	----------
	bamfile : str
		Input ailgnment file. Only Support BAMPE in current version.
	outfile : str
		Output BED file
	informat : str
		Format of the alignment file. Only Support "BAMPE" or "BAM" in current version.
	gsize : str
		Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs'
		for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7)
		and 'dm' for fruitfly (1.2e8). DEFAULT:hs
	verbose : int
		Set verbose level. 0: only show critical message, 1: show additional
		warning message, 2: show process information, 3: show debug messages.

	'''
	macs2_cmd = find_macs2()
	macs2_cmd += ' filterdup '
	dedup_cmd = "%s --keep-dup %d  --verbose %d -f %s -g %s -i %s -o %s" % (macs2_cmd, keepdup, verbose, informat, gsize, bamfile, outfile)

	if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
		logging.warning("\t\"%s\" exists and non-empy, skip deduplicating." % outfile)
	else:
		logging.info("\tRun command: %s" % dedup_cmd)
		subprocess.call(dedup_cmd, shell=True)
	fragment_count = subprocess.check_output(["wc","-l", outfile]).decode("UTF-8").split()[0]
	return (int(fragment_count))



def pileup_BED(bedfile, stype, outfile, layout, extension):
	'''
	Generate pileup track. For the ChIP sample, PE and SE data will be piled up differently.
	PE ChIP data will be piled up uisng "BEDPE". SE ChIP data will be piled up by
	extending reads to downstream direction by "extension" size. Control sampple is
	always piled up using the "SE" mode, so if the control sample is PE sequencing, it needs
	to do these conversions: BAMPE -> BEDPE -> BEDSE. SE control data will be piled
	up by extending reads to both direction by "extension/2".


	parameters
	----------

	bedfile : str
		input BED file.
	stype : str
		Sample type. Must be "ChIP" or "Ctrl"
	outfile : str
		Output bedgraph file
	layout : str
		Sequencing layout. Must be "PE" or "SE".
	extension : int
		For SE data, set as fragment size. For PE data, set as None.
	'''
	macs2_cmd = find_macs2()
	macs2_cmd += ' pileup '

	if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip pileup step." % outfile)
	else:
		if stype == 'ChIP':
			if layout == 'PE':
				pileup_cmd =  "%s -f BEDPE -i %s -o %s --verbose 0" % (macs2_cmd, bedfile, outfile)
			elif layout == 'SE':
				pileup_cmd =  "%s -f BED -i %s -o %s --extsize %d --verbose 0" % (macs2_cmd, bedfile, outfile, extension)
		elif stype == 'Ctrl':
			pileup_cmd =  "%s -f BED -B -i %s -o %s --extsize %d --verbose 0" % (macs2_cmd, bedfile, outfile, int(extension/2))
		logging.info("\tRun command: %s" %  pileup_cmd)
		subprocess.call(pileup_cmd, shell=True)


def pileup_ctrl(bedfile, out_bg_file, out_bg_norm_file, d, window_size, Ctrl_layout, sf = 1.0):
	'''
	Generate "slocal" (1kb) and "llocal" (10Kb) background tracks for the control
	sample.

	parameters
	----------

	bedfile : str
		input dedeuplicated bed file for control (*.dedup.bed).
	out_bg_file : str
		Output bedgraph file.
	out_bg_norm_file : str
		Output normalized bedgraph file.
	d : int
		Fragment size.
	window_size : int
		Window size to estimate background bias.
	chip_layout : str
		Must be "PE" or "SE". The ChIP layout determines if we need to halve the control
		track or not
	'''
	macs2_cmd = find_macs2()
	if Ctrl_layout == 'PE':
		sfactor = float(d / window_size) * 0.5 * sf
	elif Ctrl_layout == 'SE':
		sfactor = float(d / window_size) * sf
	else:
		logging.error("\tThe layout must be one of 'PE' or 'SE'.")
		sys.exit()


	if os.path.exists(out_bg_file) and os.path.getsize(out_bg_file) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip pileup step." % out_bg_file)
		pass
	else:
		pileup_cmd =  "%s pileup -f BED -B --extsize %d -i %s -o %s --verbose 0" % (macs2_cmd, int(window_size*0.5), bedfile, out_bg_file)
		logging.info("\tRun command: %s" %  pileup_cmd)
		subprocess.call(pileup_cmd, shell=True)

	if os.path.exists(out_bg_norm_file) and os.path.getsize(out_bg_norm_file) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip bdg normalization step." % out_bg_norm_file)
		pass
	else:
		normalize_cmd = "%s bdgopt -i %s -m multiply -p %f -o %s" % (macs2_cmd, out_bg_file, sfactor, out_bg_norm_file)
		logging.info("\tRun command: %s" %  normalize_cmd)
		subprocess.call(normalize_cmd, shell=True)

def max_background(bg_1K, bg_10K, bg_d, genome_noise, outfile):
	"""
	Combine and generate the maximum background noise track (bdg file) for the control sample.

	Parameters
	----------
	bg_1K : str
		The slocl background (1 Kb) bdg file.
	bg_10K : str
		The llocl background (10 Kb) bdg file.
	bg_d : str
		The d background bdg file.
	genome_noise : flaot
		Genome level background.
	outfile : str
		BEDGRAPH file containing the raw local bias from control data.
	Returns
	-------
	None.
	"""
	macs2_cmd = find_macs2()

	if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
		logging.warning("\tThe background noise track exist and non-empty: %s. Skip." % outfile)
		pass
	else:
		tmp_1K_10K = ''.join([random.choice(string.ascii_uppercase) for _ in range(16)])
		cmd_1K_vs_10K = "%s bdgcmp -m max -t %s -c %s -o %s" % (macs2_cmd, bg_1K, bg_10K, tmp_1K_10K)
		logging.info("\tRun command: %s" % cmd_1K_vs_10K)
		subprocess.call(cmd_1K_vs_10K, shell=True)

		tmp_1K_10K_d = ''.join([random.choice(string.ascii_uppercase) for _ in range(16)])
		cmd_1K_vs_10K_d = "%s bdgcmp -m max -t %s -c %s -o %s" % (macs2_cmd, tmp_1K_10K, bg_d, tmp_1K_10K_d)
		logging.info("\tRun command: %s" % cmd_1K_vs_10K_d)
		subprocess.call(cmd_1K_vs_10K_d, shell=True)


		cmd_1K_vs_10K_d_genome = "%s bdgopt -i %s -m max -p %f -o %s" % (macs2_cmd, tmp_1K_10K_d, genome_noise, outfile)
		logging.info("\tRun command: %s" % cmd_1K_vs_10K_d_genome)
		subprocess.call(cmd_1K_vs_10K_d_genome, shell=True)
	try:
		os.remove(tmp_1K_10K)
		os.remove(tmp_1K_10K_d)
	except:
		pass

def scale_bdg(in_bdg, out_bdg, sf):
	'''
	Scale bedgraph file with scaling factor (sf)

	Parameters
	----------
	in_bdg : str
		Input bedgraph file.
	out_bdg : str
		Output bedgraph file.
	sf : float
		scaling factor.

	Returns
	-------
	None.
	'''
	macs2_cmd = find_macs2()
	macs2_cmd += ' bdgopt '
	if os.path.exists(out_bdg) and os.path.getsize(out_bdg) > 0:
		logging.warning("\tThe output bedgraph file exist and non-empty:\"%s\". Skip scaling." % out_bdg)
		pass
	else:
		scale_cmd = '%s -i %s -m multiply -p %f -o %s' % (macs2_cmd, in_bdg, sf, out_bdg)
		logging.info("\tRun command: %s" % scale_cmd)
		subprocess.call(scale_cmd, shell=True)

def chip_vs_lambda(chip_file, lambda_file, out_file, method='qpois', pseudocount=0):
	'''
	Generate qvalue track

	Parameters
	----------
	chip_file : str
		Bedgraph file of ChIP sample.
	lambda_file : str
		Bedgraph file of control sample
	out_file : str
		Output bedgraph file containing qvalue at each position.
	method : str
		Method to use while calculating a score in any bin by comparing treatment value and control value. Available choices are: ppois,
		qpois, subtract, logFE, logLR, and slogLR. They represent Poisson Pvalue (-log10(pvalue) form) using control as lambda and treatment
		as observation, q-value through a BH process for poisson pvalues, subtraction from treatment, linear scale fold enrichment, log10
		fold enrichment(need to set pseudocount), log10 likelihood between ChIP-enriched model and open chromatin model(need to set
		pseudocount), symmetric log10 likelihood between two ChIP-enrichment models, or maximum value between the two tracks. Default option
		is qpois.
	Returns
	-------
	None.

	'''
	macs2_cmd = find_macs2()
	macs2_cmd += ' bdgcmp '
	if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
		logging.warning("\tThe output qvalue bedgraph file exist and non-empty:\"%s\". Skip." % out_file)
		pass
	else:
		compare_cmd = '%s -t %s -c %s  -m %s -p %d -o %s' % (macs2_cmd, chip_file, lambda_file, method, pseudocount, out_file)
		logging.info('\tRun command: %s' % compare_cmd)
		subprocess.call(compare_cmd, shell=True)

def narrow_call(in_bdg, out_bed, gap_size, min_peak_size, cut_peak = 0.05):
	'''
	Call narrow peaks.

	Parameters
	----------
	in_bdg : str
		Input qvalue bedgraph file.
	out_bed : str
		Output peak file.
	gap_size : int
		Usually set to read length. If the gap between two nearby peaks is smaller than this value, the two nearby
		peaks will be merged.
	min_peak_size : int
		Mininum peak size. Ususally set to fragment size.
	q_cut : float
		Qvalue cutoff. default = 0.05

	Returns
	-------
	None.

	'''
	macs2_cmd = find_macs2()
	macs2_cmd += ' bdgpeakcall '
	cut_peak = -log10(cut_peak)
	if os.path.exists(out_bed) and os.path.getsize(out_bed) > 0:
		logging.warning("\tThe output peak file exist and non-empty:\"%s\". Skip peak calling." % out_bed)
		pass
	else:
		callpeak_cmd = '%s --no-trackline -i %s -c %f -l %d -g %d -o %s' % (macs2_cmd, in_bdg, cut_peak, min_peak_size, gap_size, out_bed)
		logging.info('\tRun command: %s' %  callpeak_cmd)
		subprocess.call(callpeak_cmd, shell=True)


def broad_call(in_bdg, outfile, cut_peak=0.05, cut_link=0.1, min_len = 200, max_gap = 50, max_link = 800):
	'''
	Call broad peaks.

	Parameters
	----------
	in_bdg : str
		Input qvalue bedgraph file.
	outfile : float
		Output bed file.
	cut_peak : float
		Cutoff score for pvalue or qvalue in peak regions.
	cut_link : float
		Cutoff score for qvalue in linking regions.
	min_len : int
		Minimum length of peak. Usually set to d (fragmetn size). default = 200
	max_gap : int
		 Maximum gap between significant peaks. Usually set to tag size (read length). default=50
	max_link : int
		Mxximum linking between significant peaks. Usuallly set to 4*d. default = 800

	Returns
	-------
	None.

	'''
	macs2_cmd = find_macs2()
	macs2_cmd += ' bdgbroadcall '
	cut_peak = -log10(cut_peak)
	cut_link = -log10(cut_link)
	if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
		logging.warning("\tThe output peak file exist and non-empty:\"%s\". Skip peak calling." % outfile)
		pass
	else:
		callpeak_cmd = '%s --no-trackline -i %s -c %f -C %f -l %d -g %d -G %d -o %s' % (macs2_cmd, in_bdg, cut_peak, cut_link, min_len, max_gap, max_link, outfile)
		logging.info('\tRun command: %s' %  callpeak_cmd)
		subprocess.call(callpeak_cmd, shell=True)

def refine_peak(peak_bed, tag_bed, outfile):
	"""
	Refine peak to find the summit position.

	Parameters
	----------
	peak_bed : str
		Candidate peak file in BED format.
	tag_bed : str
		ChIP-seq alignment file in BED format (single end).
	outfile : str
		Output file name.

	Returns
	-------
	None.

	"""
	macs2_cmd = find_macs2()
	macs2_cmd += ' refinepeak '
	if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
		logging.warning("\tThe output peak summit file exist and non-empty:\"%s\". Skip." % outfile)
		pass
	else:
		refine_peak_cmd = '%s --verbose=0 -b %s -i %s -f BED -c 5 -w 100  -o %s' % (macs2_cmd, peak_bed, tag_bed, outfile)
		logging.info('\tRun command: %s' %  refine_peak_cmd)
		subprocess.call(refine_peak_cmd, shell=True)

def estimate_d(bed_file, r_file, informat = 'BED', gsize='hs', mfold='5  50'):
	'''
	Parameters
	----------
	bed_file : str
		Input BED fil. MACS2 'predictd' function only works for single end data.

	Returns
	-------
	None.

	'''
	estimated_d = 0
	macs2_cmd = find_macs2()
	macs2_cmd += ' predictd '
	rfile = r_file
	predict_d_cmd = '%s --verbose=0 -i %s -f %s  -g %s  -m %s --rfile %s' % (macs2_cmd, bed_file, informat, gsize, mfold, rfile)
	logging.info('\tRun command: %s' %  predict_d_cmd)
	subprocess.call(predict_d_cmd, shell=True)
	#legend('right','alt lag(s) : 137',bty='n')

	if os.path.exists(rfile) and os.path.getsize(rfile) > 0:
		for l in open(rfile,'r'):
			l = l.strip()
			if not l.startswith("legend('right'"):
				continue
			f = l.split(',')
			#print (f[1])
			estimated_d = f[1].replace("'","").replace("alt lag(s) : ", "")
	else:
		logging.warning('\tCannot estimate fragment size. Set to 146')
		estimated_d = 146

	return int(estimated_d)

