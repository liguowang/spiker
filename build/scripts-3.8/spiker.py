#!python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 11:02:30 2019

@author: m102324
"""

import sys
import logging
import string
import random
import os
from optparse import OptionParser
from spikers import run_bowtie2,run_macs2,process_bam,PE,write_bw


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.3"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():

	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)

	parser.add_option("--t1",action="store",type="string",dest="chip_R1",help="FASTQ file (read1) for ChIP sample. Can be regular plain text file or compressed file (.gz, .bz2). Mutually exclusive with '-t'.")
	parser.add_option("--t2",action="store",type="string",dest="chip_R2",help="FASTQ file (reas2) for ChIP sample. Can be regular plain text file or compressed file (.gz, .bz2). Mutually exclusive with '-t'. Ignore this for single-end sequencing.")
	parser.add_option("-t", "--treat",action="store",type="string",dest="chip_bam",help="BAM file of ChIP sample. The BAM file must be sorted and indexed. Mutually exclusive with '--t1' and '--t2'. ")
	parser.add_option("--c1",action="store",type="string",dest="ctrl_R1",help="FASTQ file (read1) for Control sample. Can be regular plain text file or compressed file (.gz, .bz2). Mutually exclusive with '-c'.")
	parser.add_option("--c2",action="store",type="string",dest="ctrl_R2",help="FASTQ file (reas2) for Control sample. Can be regular plain text file or compressed file (.gz, .bz2). Mutually exclusive with '-c'. Ignore this for single-end sequencing.")
	parser.add_option("-c", "--control", action="store",type="string",dest="ctrl_bam",help="BAM file of Control sample. Mutually exclusive with '--c1' and '--c2'. The BAM file must be sorted and indexed.")
	parser.add_option("-o","--output", action="store",type="string",dest="outfile",help="Prefix of output files.")
	parser.add_option("--bt2-index",action="store",type="string",dest="bt2_index",help="The prefix (minus trailing .X.bt2) for bowtie2 index files. Ignore this option if BAM files were provided by '-t' and '-c'.")
	parser.add_option("-n",action="store",type='int', default=1000000, dest="n_reads",help="Number of alignments from the BAM file used to tell the sequencing layout (PE or SE), and estiamte the fragment size 'd'. default=%default")
	parser.add_option("-g","--genome-size",action="store",type='string', default='hs', dest="g_size",help="Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8). default=%default")
	parser.add_option("-p","--proc",action="store",type='int', default=8, dest="n_threads",help="Number of threads. default=%default")
	parser.add_option("--mfold",action="store",dest="m_fold", default='5 50', help="Select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. Fold-enrichment in regions must be lower than upper limit, and higher than the lower limit. Use as \"-m 10 30\". DEFAULT:5 50")
	parser.add_option("--spikeIn",action="store_true",dest="spike_in", default=False, help="Set this flag if ChIP and control samples contains exogenous reads as splike-in. Please note, you also need to specify --tsf and --csf.")
	parser.add_option("--tsf",action="store",type='float', dest="treat_sf", help="Scaling factor for treatment. This will be applied to the pileup bedgraph file of treatment (*.treat.pileup.bdg).")
	parser.add_option("--csf",action="store",type='float', dest="control_sf", help="Scaling factor for control. This will be applied to the pileup bedgraph file of maximum background (*.control.pileup.max.bdg).")
	parser.add_option("--q-peak",action="store",type='float', dest="q_cutoff", default=0.05, help="Qvalue cutoff for peaks. default=%default")
	parser.add_option("--q-link",action="store",type='float', dest="q_link_cut", default=0.1, help="Qvalue cutoff for linking regions. default=%default")
	parser.add_option("--bw",action="store_true",dest="bigwig", default=False, help="If set, generate bigwig files for ChIP pileup and control pileup.")
	parser.add_option("--maxgap",action="store", type='int', dest="max_gap", default=100, help="maximum gap between significant points in a peak. default=%default")
	parser.add_option("--broad",action="store_true",dest="broad_peak", default=False, help="If set, call broad peaks.")
	parser.add_option("--frip",action="store_true",dest="frip", default=False, help="If set, calculate FRiP (the Fraction of Reads In called Peaks) score using the BAM and peak files.")
	parser.add_option("--cleanup",action="store_true",dest="clean_up", default=False, help="If set, clean up the intermediate files. When not set, intermediate files are kept so that rerun the workflwo will be much faster.")
	parser.add_option("--refine",action="store_true",dest="refine_peak", default=False, help="If set, detect peak summit position.")
	parser.add_option("--verbose",action="store_true",dest="debug", default=False, help="If set, print detailed information for debugging.")

	(options,args)=parser.parse_args()

	#DEGUB->INFO->WARNING->ERROR->CRITICAL
	if options.debug:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
	else:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)


	#=========================================================================
	# Set up genome size
	#=========================================================================
	genome_size = 0
	gsizes = {'hs':2.7e9, 'mm':1.87e9, 'ce':9e7, 'dm':1.2e8}
	if options.g_size in gsizes:
		genome_size = gsizes[options.g_size]
	else:
		genome_size = int(options.g_size)
	if genome_size <= 0:
		logging.error("Unknown genome size: %d" % genome_size)
		sys.exit()

	#=========================================================================
	# Check spike-in options
	#=========================================================================
	if options.spike_in:
		if not (options.treat_sf):
			parser.print_help()
			print('\nError:\n--tsf must be specified for spike-in data.', file=sys.stderr)
			sys.exit()
		if not (options.control_sf):
			parser.print_help()
			print('\nError:\n--csf must be specified for spike-in data.', file=sys.stderr)
			sys.exit()
	if options.treat_sf and options.control_sf:
		options.spike_in = True


	logging.info("Running ChIP-seq workflow ...")


	#=========================================================================
	# Input single end fastq files
	#=========================================================================
	if all([options.chip_R1, options.ctrl_R1, options.outfile, options.bt2_index])  and not any([options.chip_R2, options.ctrl_R2, options.chip_bam, options.ctrl_bam]):
		bam_input = False
		for test_file in (options.chip_R1, options.ctrl_R1):
			if (os.path.exists(test_file) and os.path.getsize(test_file) > 0):
				pass
			else:
				parser.print_help()
				print('\nError:\nCannot find "%s".' % test_file, file=sys.stderr)
				sys.exit()
		logging.info("Input single-end FASTQ files ...\n")
		logging.info("Step_1: Map reads to the reference genome")
		#bowtie2 for ChIP
		logging.info("Step_1.1: Map ChIP SE reads to the reference genome with bowtie2")
		chip_bam = options.outfile + '.treat.bam'
		chip_bam_sorted = options.outfile + '.treat.sorted.bam'
		run_bowtie2.bowtie2_map(bt2_prefix = options.bt2_index, fq1 = options.chip_R1, fq2 = None, outbam = chip_bam, nthread = options.n_threads)

		#bowtie2 for control
		logging.info("Step1.2: Map control SE reads to the reference genome with bowtie2")
		ctrl_bam = options.outfile + '.control.bam'
		ctrl_bam_sorted = options.outfile + '.control.sorted.bam'
		run_bowtie2.bowtie2_map(bt2_prefix = options.bt2_index, fq1 = options.ctrl_R1, fq2 = None, outbam = ctrl_bam, nthread = options.n_threads)


	#=========================================================================
	# Input paired end fastq files
	#=========================================================================
	elif all([options.chip_R1, options.ctrl_R1, options.chip_R2, options.ctrl_R2, options.outfile, options.bt2_index])  and not any([options.chip_bam, options.ctrl_bam is None]):
		bam_input = False
		for test_file in (options.chip_R1, options.chip_R2, options.ctrl_R1, options.ctrl_R):
			if (os.path.exists(test_file) and os.path.getsize(test_file) > 0):
				pass
			else:
				parser.print_help()
				print('\nError:\nCannot find "%s".' % test_file, file=sys.stderr)
				sys.exit()
		logging.info("Input pair-end FASTQ files ...\n")
		logging.info("Step_1: Map reads to the reference genome")
		#bowtie2 for ChIP
		logging.info("Step_1.1: Map ChIP PE reads to the reference genome with bowtie2")
		chip_bam = options.outfile + '.treat.bam'
		chip_bam_sorted = options.outfile + '.treat.sorted.bam'
		run_bowtie2.bowtie2_map(bt2_prefix = options.bt2_index, fq1 = options.chip_R1, fq2 = options.chip_R2, outbam = chip_bam, nthread = options.n_threads)

		#bowtie2 for control
		logging.info("Step1.2: Map control PE reads to the reference genome with bowtie2")
		ctrl_bam = options.outfile + '.control.bam'
		ctrl_bam_sorted = options.outfile + '.control.sorted.bam'
		run_bowtie2.bowtie2_map(bt2_prefix = options.bt2_index, fq1 = options.ctrl_R1, fq2 = options.ctrl_R2, outbam = ctrl_bam, nthread = options.n_threads)


	#=========================================================================
	# Input BAM files
	#=========================================================================
	elif not any([options.chip_R1, options.ctrl_R1, options.chip_R2, options.ctrl_R2]) and all([options.outfile, options.chip_bam, options.ctrl_bam]):
		bam_input = True
		for test_file in (options.chip_bam, options.ctrl_bam, options.chip_bam + '.bai', options.ctrl_bam + '.bai'):
			if (os.path.exists(test_file) and os.path.getsize(test_file) > 0):
				pass
			else:
				parser.print_help()
				print('\nError:\nCannot find "%s".' % test_file, file=sys.stderr)
				sys.exit()
		logging.info("Input BAM files. Skip Step_1 (reads mapping).")
		chip_bam_sorted = options.chip_bam
		ctrl_bam_sorted = options.ctrl_bam
	else:
		parser.print_help()
		print ("\nUsage examples:\n", file=sys.stderr)
		print ("\n#Input single-end fastq files:", file=sys.stderr)
		print ("\t$ spiker.py --t1 H3K27ac.fastq.gz --c1 control.fastq.gz --bt2-index /data/GRCh38 --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac", file=sys.stderr)
		print ("\t$ spiker.py --broad --t1 H3K27ac.fastq.gz --c1 control.fastq.gz --bt2-index /data/GRCh38 --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac", file=sys.stderr)

		print ("\n#Input paired-end fastq files:", file=sys.stderr)
		print ("\t$ spiker.py --t1 H3K27ac_R1.fastq.gz --t2 H3K27ac_R2.fastq.gz --c1 control_R1.fastq.gz --c2 control_R2.fastq.gz --bt2-index /data/GRCh38 --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac", file=sys.stderr)
		print ("\t$ spiker.py --broad --t1 H3K27ac_R1.fastq.gz --t2 H3K27ac_R2.fastq.gz --c1 control_R1.fastq.gz --c2 control_R2.fastq.gz --bt2-index /data/GRCh38 --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac", file=sys.stderr)

		print ("\n#Input BAM files:", file=sys.stderr)
		print ("\t$ spiker.py -t H3K27ac.sorted.bam -c control.sorted.bam --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac", file=sys.stderr)
		print ("\t$ spiker.py --broad -t H3K27ac.sorted.bam -c control.sorted.bam --spikeIn --csf 1.23 --tsf 0.95 -o H3K27ac", file=sys.stderr)
		print ()
		sys.exit()


	#=========================================================================
	# Extract some information from ChIP and control BAM files, include:
	#   * layout (PE or SE)
	#   * read length
	#   * total mapped reads
	#   * chrom sizes
	#   * fragment size (only for BAM-PE). fragment size cannot be calculated
	#     from the BAM-SE.
	#=========================================================================
	logging.info("Step_2: Extract information from BAM files ...")

	#check sequencing layout
	logging.info("Step_2.1: Check sequencing layout from BAM file ...")
	chip_layout = process_bam.get_layout(bamfile = chip_bam_sorted, n=options.n_reads)
	logging.info("\tThe layout of %s is %s" % (chip_bam_sorted, chip_layout ))
	ctrl_layout = process_bam.get_layout(bamfile = ctrl_bam_sorted, n=options.n_reads)
	logging.info("\tThe layout of %s is %s" % (ctrl_bam_sorted, ctrl_layout ))
	if chip_layout == 'Mix' or chip_layout == 'Unknown':
		logging.error("\tInput ChIP data must be PE or SE. Exit")
		sys.exit()
	if ctrl_layout == 'Mix' or ctrl_layout == 'Unknown':
		logging.error("\tInput control data must be PE or SE. Exit")
		sys.exit()

	# calculate fragment size (d), read length, chrom sizes
	logging.info("Step_2.2: Extract information from ChIP BAM files ...")
	d_chip, chip_readSize, chrom_sizes = process_bam.bam_info(chip_bam_sorted, n = options.n_reads, layout = chip_layout)
	if d_chip:
		logging.info("\tThe average fragment size of ChIP sample is: %d" % int(d_chip))
	chip_total_mapped = process_bam.total_mapped(chip_bam_sorted)
	logging.info("\tThe total mapped reads in ChIP sample is: %d" % chip_total_mapped)

	logging.info("Step_2.3: Extract information from control BAM files ...")
	d_ctrl, ctrl_readSize, chrom_sizes = process_bam.bam_info(ctrl_bam_sorted, n = options.n_reads, layout = ctrl_layout)
	if d_ctrl:
		logging.info("\tThe average fragment size of control sample is: %d" % int(d_ctrl))
	ctrl_total_mapped = process_bam.total_mapped(ctrl_bam_sorted)
	logging.info("\tThe total mapped reads in control sample is: %d" % ctrl_total_mapped)



	#=========================================================================
	# Deduplication and convert BAM into BED file (for ChIP sample)
	#=========================================================================
	logging.info("Step_3: Deduplication ...")
	logging.info("Step_3.1: Remove duplidate reads from ChIP BAM file ...")
	chip_dedup_bed_PE =  options.outfile + '.treat.dedup.PE.bed'
	chip_dedup_bed_SE =  options.outfile + '.treat.dedup.SE.bed'

	if chip_layout == 'PE':
		chip_frag_num = run_macs2.deduplicate(bamfile = chip_bam_sorted, outfile = chip_dedup_bed_PE, informat = 'BAMPE', gsize = options.g_size)
	elif chip_layout == 'SE':
		chip_predictd_r = options.outfile + '.treat.predict_d.r'
		chip_frag_num = run_macs2.deduplicate(bamfile = chip_bam_sorted, outfile = chip_dedup_bed_SE, informat = 'BAM', gsize = options.g_size)
		# for SE data, have to re-estimate fragmetn size
		d_chip = run_macs2.estimate_d(bed_file = chip_dedup_bed_SE, r_file = chip_predictd_r, informat = 'BED', gsize = options.g_size, mfold = options.m_fold)
		logging.info("\tThe average fragment size of ChIP sample is: %d" % int(d_chip))

	logging.info("\tNumber of unique ChIP fragments: %d" % chip_frag_num)
	if options.spike_in and options.treat_sf != 1:
		chip_frag_num = int(chip_frag_num * options.treat_sf)
		logging.info("\tNumber of unique ChIP fragments after spike-in adjustment: %d" % chip_frag_num)


	#=========================================================================
	# Deduplication and convert BAM into BED file (for control sample)
	#=========================================================================
	logging.info("Step_3.2: Remove duplidate reads from control BAM file ...")
	ctrl_dedup_bed_PE =  options.outfile + '.control.dedup.PE.bed'
	ctrl_dedup_bed_SE =  options.outfile + '.control.dedup.SE.bed'
	if ctrl_layout == 'PE':
		ctrl_frag_num = run_macs2.deduplicate(bamfile = ctrl_bam_sorted, outfile = ctrl_dedup_bed_PE, informat = 'BAMPE', gsize = options.g_size)
	elif ctrl_layout == 'SE':
		ctrl_predictd_r = options.outfile + '.control.predict_d.r'
		ctrl_frag_num = run_macs2.deduplicate(bamfile = ctrl_bam_sorted, outfile = ctrl_dedup_bed_SE, informat = 'BAM', gsize = options.g_size)
		# for SE data, have to re-estimate fragmetn size
		d_ctrl = run_macs2.estimate_d(bed_file = ctrl_dedup_bed_SE, r_file = ctrl_predictd_r, informat = 'BED', gsize = options.g_size, mfold = options.m_fold)
		logging.info("\tThe average fragment size of control sample is: %d" % int(d_ctrl))

	logging.info("\tNumber of unique control fragments: %d" % ctrl_frag_num)
	if options.spike_in and options.control_sf != 1:
		ctrl_frag_num = int(ctrl_frag_num * options.control_sf)
		logging.info("\tNumber of unique control fragments after spike-in adjustment: %d" % ctrl_frag_num)


	#=========================================================================
	# If input is PE data, convert deduplicated BED-PE into BED-SE
	# ChIP BED-SE file is to "refine peaks" (find summit)
	# Control BED-SE file is to build d-pileup, slocal-pileup, llocal-pileup
	# Nothing will be done if the input is SE data
	#=========================================================================
	logging.info("Step_4: Convert PE to SE")
	if chip_layout == 'PE':
		logging.info("Step_4.1: Convert BED (PE) into BED (SE) for the ChIP sample ...")
		PE.pe2se(in_PE = chip_dedup_bed_PE, out_SE = chip_dedup_bed_SE, rl = chip_readSize)
	else:
		logging.info("Step_4.1: Convert BED (PE) into BED (SE) for the ChIP sample (Skip)")
		pass
	if ctrl_layout == 'PE':
		logging.info("Step_4.2: Convert BED (PE) into BED (SE) for the control sample ...")
		PE.pe2se(in_PE = ctrl_dedup_bed_PE, out_SE = ctrl_dedup_bed_SE, rl = ctrl_readSize)
	else:
		logging.info("Step_4.2: Convert BED (PE) into BED (SE) for the control sample (Skip)")
		pass


	#=========================================================================
	# Create the pileup track file (in bedGraph format) for ChIP sample
	# For PE data, piled up using BED-PE mode
	# for SE data, extend 'd' to downstream
	# The resulting bedGraph file will be scaled if --tsf was specified
	# The resulting bedGraph file will be used for peak calling
	#=========================================================================
	logging.info("Step_5: Pileup ChIP alignments ...")
	chip_pileup =  options.outfile + '.treat.pileup.bdg'
	#run_macs2.pileup_BED_PE(bedpefile = chip_dedup_bed, outfile = chip_pileup)
	if chip_layout == 'PE':
		run_macs2.pileup_BED(bedfile = chip_dedup_bed_PE, stype = "ChIP", outfile = chip_pileup, layout = 'PE', extension = None)
	elif chip_layout == 'SE':
		run_macs2.pileup_BED(bedfile = chip_dedup_bed_SE, stype = "ChIP", outfile = chip_pileup, layout = 'SE', extension = d_chip)
	if options.spike_in and options.treat_sf != 1:
		logging.info("Step_5.1: Scale ChIP pileup uisng Spike-In ...")
		chip_pileup_spikeIn = options.outfile + '.treat.pileup.SpikeIn_scaled.bdg'
		run_macs2.scale_bdg(in_bdg = chip_pileup, out_bdg = chip_pileup_spikeIn, sf = options.treat_sf)

	#=========================================================================
	# Create the pileup track file (in bedGraph format) for control sample
	# PE data will be converted into SE
	# for SE data, extend 'd/2' to both downstream and upstream
	# This is called d-pileup
	# The resulting bedGraph file will be halved if control data is PE
	#=========================================================================

	logging.info("Step_6: Pileup control alignments (d-pileup) ...")
	ctrl_d_pileup =  options.outfile + '.control.pileup.d.bdg'
	run_macs2.pileup_BED(bedfile = ctrl_dedup_bed_SE,  stype = "Ctrl", outfile = ctrl_d_pileup, layout = 'SE', extension = d_ctrl)

	if ctrl_layout == 'PE':
		ctrl_d_pileup_spikeIn = options.outfile + '.control.pileup.d.SpikeIn_scaled.bdg'
		if options.spike_in and options.control_sf != 1:
			logging.info("Step_6.1: Halve control's d-pileup and scale it uisng Spike-In ...")
			run_macs2.scale_bdg(in_bdg = ctrl_d_pileup, out_bdg = ctrl_d_pileup_spikeIn, sf = options.control_sf * 0.5)
		else:
			logging.info("Step_6.1: Halve control's d-pileup ...")
			tmp_file = ''.join([random.choice(string.ascii_uppercase) for _ in range(10)])
			run_macs2.scale_bdg(in_bdg = ctrl_d_pileup, out_bdg = tmp_file, sf = 0.5)
			os.rename(tmp_file, ctrl_d_pileup)
	else:
		ctrl_d_pileup_spikeIn = options.outfile + '.control.pileup.d.SpikeIn_scaled.bdg'
		if options.spike_in and options.control_sf != 1:
			logging.info("Step_6.1: Scale control uisng Spike-In ...")
			run_macs2.scale_bdg(in_bdg = ctrl_d_pileup, out_bdg = ctrl_d_pileup_spikeIn, sf = options.control_sf)


	#=========================================================================
	# Create the slocal track file for control sample
	# PE data will be converted into SE
	# for SE data, extend 'window/2' to both downstream and upstream (default, window = 1 Kb)
	# This is called slocal-pileup
	# The resulting bedGraph file will be halved if control data is PE
	# The resulting bedGraph file will be multiply factor provided by --csf
	# The resulting bedGraph file will be multiply factor d/window
	#=========================================================================

	logging.info("Step_6.2: Build slocal-pileup ...")
	ctrl_1Kb_pileup =  options.outfile + '.control.pileup.1Kb.bdg'
	ctrl_1Kb_norm_pileup =  options.outfile + '.control.pileup.1Kb_norm.bdg'
	run_macs2.pileup_ctrl(bedfile = ctrl_dedup_bed_SE, out_bg_file = ctrl_1Kb_pileup, out_bg_norm_file = ctrl_1Kb_norm_pileup, d = d_ctrl, window_size = 1000, Ctrl_layout = ctrl_layout, sf = options.control_sf)

	#=========================================================================
	# Create the llocal track file for control sample
	# PE data will be converted into SE
	# for SE data, extend 'window/2' to both downstream and upstream (default, window = 10 Kb)
	# This is called llocal-pileup
	# The resulting bedGraph file will be halved if control data is PE
	# The resulting bedGraph file will be multiply factor provided by --csf
	# The resulting bedGraph file will be multiply factor d/window
	#=========================================================================
	logging.info("Step_6.3: Build llocal-pileup ...")
	ctrl_10Kb_pileup =  options.outfile + '.control.pileup.10Kb.bdg'
	ctrl_10Kb_norm_pileup =  options.outfile + '.control.pileup.10Kb_norm.bdg'
	run_macs2.pileup_ctrl(bedfile = ctrl_dedup_bed_SE, out_bg_file = ctrl_10Kb_pileup, out_bg_norm_file = ctrl_10Kb_norm_pileup, d = d_ctrl, window_size = 10000, Ctrl_layout = ctrl_layout, sf = options.control_sf)


	#=========================================================================
	# Take the maximum out of (d-pileup, slocal-pileup, llocal-pileup, and
	# genome-background)
	# ctrl_maxNoise_pileup is used for peak calling
	#=========================================================================

	# genome level background
	genome_background = (ctrl_frag_num * d_ctrl)/genome_size
	logging.info("Step_6.4: Estimated genome background: %f" % genome_background)
	#generate background noise track
	logging.info("Step_6.5: Generating the maximum background noise track from (d-pileup, slocal-pileup, llocal-pileup and genome_background) ...")
	ctrl_maxNoise_pileup = options.outfile + '.control.pileup.max.bdg'
	run_macs2.max_background(bg_1K = ctrl_1Kb_norm_pileup, bg_10K = ctrl_10Kb_norm_pileup, bg_d = ctrl_d_pileup_spikeIn, genome_noise = genome_background, outfile = ctrl_maxNoise_pileup)

	#=========================================================================
	# Peak calling (spike-in ChIP-seq)
	#=========================================================================
	if options.spike_in:
		#compare ChIP to lambda
		logging.info("Step_7: Compare ChIP to control to generate qvalue track (use qpois) ...")
		chip_qvalue_bdg = options.outfile + '.qvalue.bdg'
		run_macs2.chip_vs_lambda(chip_file = chip_pileup_spikeIn, lambda_file = ctrl_maxNoise_pileup, out_file = chip_qvalue_bdg)

		# call peak
		if options.broad_peak:
			logging.info("Step_8: Calling broad peaks ...")
			broad_peak = options.outfile + '.broadPeak'
			run_macs2.broad_call(in_bdg = chip_qvalue_bdg, outfile = broad_peak, cut_peak= options.q_cutoff, cut_link = options.q_link_cut, min_len = d_chip, max_gap = max(chip_readSize, options.max_gap), max_link = d_chip*4)
			if options.refine_peak:
				broad_peak_summit = options.outfile + '.broadPeak.summit.bed'
				logging.info("Step_8.1: refine peaks ...")
				run_macs2.refine_peak(peak_bed = broad_peak, tag_bed = chip_dedup_bed_SE, outfile = broad_peak_summit)
		else:
			logging.info("Step_8: Calling narrow peaks ...")
			narrow_peak = options.outfile + '.narrowPeak'
			run_macs2.narrow_call(in_bdg = chip_qvalue_bdg, out_bed = narrow_peak, gap_size = max(chip_readSize, options.max_gap), min_peak_size = d_chip, cut_peak = options.q_cutoff)
			if options.refine_peak:
				narrow_peak_summit = options.outfile + '.narrowPeak.summit.bed'
				logging.info("Step_8.1: refine peaks ...")
				run_macs2.refine_peak(peak_bed = narrow_peak, tag_bed = chip_dedup_bed_SE, outfile = narrow_peak_summit)
	#=========================================================================
	# Peak calling (regular ChIP-seq)
	#=========================================================================
	else:
		logging.info("Step_7: Compare ChIP to control to generate qvalue track (use qpois) ...")
		logging.info("\tScale ChIP & control to the same sequencing depth")
		chip_ctrl_ratio = chip_frag_num / ctrl_frag_num
		#no need to scale ChIP
		chip_depth_scaled = chip_pileup
		control_depth_scaled = options.outfile + '.control.pileup.depth_scaled.bdg'
		run_macs2.scale_bdg(in_bdg = ctrl_maxNoise_pileup, out_bdg = control_depth_scaled, sf = chip_ctrl_ratio)

		#compare ChIP to lambda
		chip_qvalue_bdg = options.outfile + '.qvalue.bdg'
		run_macs2.chip_vs_lambda(chip_file = chip_depth_scaled, lambda_file = control_depth_scaled, out_file = chip_qvalue_bdg )

		# call peak
		if options.broad_peak:
			logging.info("Step_8: Calling broad peaks ...")
			broad_peak = options.outfile + '.broadPeak'
			run_macs2.broad_call(in_bdg = chip_qvalue_bdg, outfile = broad_peak, cut_peak = options.q_cutoff, cut_link = options.q_link_cut, min_len = d_chip, max_gap = max(chip_readSize, options.max_gap), max_link = d_chip*4)
			if options.refine_peak:
				broad_peak_summit = options.outfile + '.broadPeak.summit.bed'
				logging.info("Step_8.1: refine peaks ...")
				run_macs2.refine_peak(peak_bed = broad_peak, tag_bed = chip_dedup_bed_SE, outfile = broad_peak_summit)
		else:
			logging.info("Step_8: Calling narrow peaks ...")
			narrow_peak = options.outfile + '.narrowPeak'
			run_macs2.narrow_call(in_bdg = chip_qvalue_bdg, out_bed = narrow_peak, gap_size = max(chip_readSize, options.max_gap), min_peak_size = d_chip, cut_peak = options.q_cutoff)
			if options.refine_peak:
				narrow_peak_summit = options.outfile + '.narrowPeak.summit.bed'
				logging.info("Step_8.1: refine peaks ...")
				run_macs2.refine_peak(peak_bed = narrow_peak, tag_bed = chip_dedup_bed_SE, outfile = narrow_peak_summit)

	#=========================================================================
	# Generate bigWig files
	#=========================================================================
	if options.bigwig:
		logging.info("Step_9: Convert bedGraph into bigWig files ...")
		if options.spike_in:
			logging.info("Step_9.1: Convert '%s' into bigwig ..." % chip_pileup_spikeIn)
			chip_pileup_spikeIn_bw = options.outfile + '.treat.pileup.SpikeIn_scaled.bigWig'
			write_bw.bdg2bw(in_bdg = chip_pileup_spikeIn, out_bw = chip_pileup_spikeIn_bw, chromSizes = chrom_sizes)

			logging.info("Step9.2: Convert '%s' into bigwig ..." % ctrl_maxNoise_pileup)
			ctrl_maxNoise_pileup_bw = options.outfile + '.control.pileup.max.bigWig'
			write_bw.bdg2bw(in_bdg = ctrl_maxNoise_pileup, out_bw = ctrl_maxNoise_pileup_bw, chromSizes = chrom_sizes)
		else:
			logging.info("Step_9.1: Convert '%s' into bigwig ..." % chip_depth_scaled)
			chip_pileup_bw = options.outfile + '.treat.pileup.bigWig'
			write_bw.bdg2bw(in_bdg = chip_depth_scaled, out_bw = chip_pileup_bw, chromSizes = chrom_sizes)

			logging.info("Step9.2: Convert '%s' into bigwig ..." % control_depth_scaled)
			ctrl_maxNoise_pileup_bw = options.outfile + '.control.pileup.bigWig'
			write_bw.bdg2bw(in_bdg = control_depth_scaled, out_bw = ctrl_maxNoise_pileup_bw, chromSizes = chrom_sizes)

	#=========================================================================
	# Calculate FRiP score
	#=========================================================================
	if options.frip:
		from hischiplib import cal_FRiP
		if options.broad_peak:
			logging.info("Step_10: Calculate FRiP score using '%s' and '%s'" % (broad_peak, chip_bam_sorted))
			frip_score = cal_FRiP.frip_score(peak_file = broad_peak, bam_file = chip_bam_sorted, nthreads = options.n_threads)
		else:
			logging.info("Step_10: Calculate FRiP score using '%s' and '%s'" % (narrow_peak, chip_bam_sorted))
			frip_score = cal_FRiP.frip_score(peak_file = narrow_peak, bam_file = chip_bam_sorted, nthreads = options.n_threads)
		logging.info("\tFRiP = %f" % frip_score)

	#=========================================================================
	# Cleaning up intermediate files
	#=========================================================================
	if options.clean_up:
		logging.info("Step_11: cleaning up intermediate files ...")
		if bam_input is False:
			tmp_files = [chip_bam, ctrl_bam, chip_dedup_bed_PE, chip_dedup_bed_SE, ctrl_dedup_bed_PE, ctrl_dedup_bed_SE, ctrl_d_pileup, ctrl_d_pileup_spikeIn, ctrl_1Kb_pileup, ctrl_1Kb_norm_pileup, ctrl_10Kb_pileup, ctrl_10Kb_norm_pileup, chip_qvalue_bdg, 'tmp_1K_10K.bdg', 'tmp_1K_10K_d.bdg']
		else:
			tmp_files = [chip_dedup_bed_PE, chip_dedup_bed_SE, ctrl_dedup_bed_PE, ctrl_dedup_bed_SE, ctrl_d_pileup, ctrl_d_pileup_spikeIn, ctrl_1Kb_pileup, ctrl_1Kb_norm_pileup, ctrl_10Kb_pileup, ctrl_10Kb_norm_pileup, chip_qvalue_bdg, 'tmp_1K_10K.bdg', 'tmp_1K_10K_d.bdg']
		for tmp in tmp_files:
			try:
				logging.debug("\t:removing %s" % tmp)
				os.remove(tmp)
			except:
				pass



if __name__=='__main__':
	main()
