#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:00:02 2020

@author: m102324
"""

import sys,os
import pysam
from optparse import OptionParser
import logging
import shutil
import subprocess


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.3"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def create_headers(bamfile, prog_name='spiker', prog_ver = __version__, co=[]):
	"""
	Create BAM headers for human and exogenous BAM files. Note, to distinguish
	the chromosomes of human and exogenous species, the chrom ID of
	exogenous chromosomes cannot start with "chr" (in other words, "chr" is
	reserved only for human chromosomes). For example, prefix "dm6_" was added
	to the chromosome IDs of Drosophila melanogaster (fruit fly).
	Parameters
	----------
	bamfile : AlignmentFile object
		a=pysam.AlignmentFile('33_K27Ac-112.sorted.bam','rb')
	Returns
	-------
	None.
	"""
	bam_header = bamfile.header #b am header from the composite genome.
	hs_header = {} # bam header for human genome alignments
	ex_header = {} # bam header for exogenous genome (such as fly) alignments
	# update 'SQ'
	for key in bam_header.keys():
		if key == 'SQ':
			if 'SQ' not in hs_header:
				hs_header['SQ'] = []
			if 'SQ' not in ex_header:
				ex_header['SQ'] = []
			lst = bam_header['SQ']
			for i in lst:
				if i['SN'].startswith('chr'):
					hs_header['SQ'].append(i)
				else:
					ex_header['SQ'].append(i)
		else:
			if key not in hs_header:
				hs_header[key] = bam_header[key]
	# update 'PG'
	if 'PG' in hs_header:
		hs_header['PG'].append( {'ID':prog_name,'PN':prog_name, 'VN':prog_ver})
	else:
		hs_header['PG'] = [{'ID':prog_name,'PN':prog_name, 'VN':prog_ver}]
	if 'PG' in ex_header:
		ex_header['PG'].append( {'ID':prog_name,'PN':prog_name, 'VN':prog_ver})
	else:
		ex_header['PG'] = [{'ID':prog_name,'PN':prog_name, 'VN':prog_ver}]
	# update 'CO'
	for comment in co:
		if 'CO' in hs_header:
			hs_header['CO'].append(comment)
		else:
			hs_header['CO'] = [comment]
		if 'CO' in ex_header:
			ex_header['CO'].append(comment)
		else:
			ex_header['CO'] = [comment]
	return(hs_header, ex_header)

def lookup_refid(hd,chr_name):
	lst = hd['SQ']
	for i in range(len(lst)):
		if chr_name == lst[i]['SN']:
			return i
	return None

def sort_bam(bamfile, threads = 1):

	sorted_bam = bamfile.replace('.bam','.sorted.bam')
	sorted_bam_bai = sorted_bam + '.bai'

	# find samtools command
	samtools_cmd = shutil.which("samtools")
	if samtools_cmd is None:
		logging.error("\tCannot find the \"samtools\" command!")
		sys.exit()
	logging.debug("\tFound samtools command:\"%s\"" %  samtools_cmd)


	# sort BAM file
	if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip BAM sorting." % sorted_bam)
		pass
	else:
		logging.info("\tSorting %s ..." % bamfile)
		samtools_sort = "%s sort -@ %d %s > %s" % (samtools_cmd, threads, bamfile, sorted_bam)
		subprocess.call(samtools_sort, shell=True)

	# index BAM
	if os.path.exists(sorted_bam_bai) and os.path.getsize(sorted_bam_bai) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip BAM indexing." % sorted_bam_bai)
		pass
	else:
		logging.info("\Indexing %s ..." % sorted_bam)
		samtools_index = "%s index %s" % (samtools_cmd, sorted_bam)
		subprocess.call(samtools_index, shell=True)


def divided_bam(bam_file, outfile, q_cut=30, chr_prefix='dm6_', threads = 1):
	unmapped_reads = 0
	qcfail_reads = 0
	duplicate_reads = 0
	secondary_reads = 0
	low_maq = 0
	diff_genome = 0	#read pairs mapped to fly and human simutaneously
	fly_reads = 0
	human_reads = 0
	samfile = pysam.AlignmentFile(bam_file,'rb')
	human_header, ex_header = create_headers(samfile)
	OUT = open(outfile + '.report.txt','w')
	HU = pysam.AlignmentFile(outfile + '_human.bam', "wb", header=human_header)
	EX = pysam.AlignmentFile(outfile + '_exogenous.bam', "wb", header=ex_header)
	BO = pysam.AlignmentFile(outfile + '_both.bam', "wb", template=samfile)
	NE = pysam.AlignmentFile(outfile + '_neither.bam', "wb", template=samfile) #unmapped, qc failed, duplicate
	logging.info("Read the BAM file: %s" % bam_file)
	#PE = False
	try:
		while(1):
			aligned_read = next(samfile)
			if aligned_read.is_unmapped:
				unmapped_reads += 1
				NE.write(aligned_read)
				continue
			elif aligned_read.is_qcfail:
				qcfail_reads += 1
				NE.write(aligned_read)
				continue
			elif aligned_read.is_duplicate:
				duplicate_reads += 1
				NE.write(aligned_read)
				continue
			elif aligned_read.is_secondary:
				secondary_reads += 1
				NE.write(aligned_read)
				continue
			elif aligned_read.mapq < q_cut:
				NE.write(aligned_read)
				low_maq += 1
				continue
			else:
				if aligned_read.is_paired:
					#PE = True
					read1_chr = samfile.get_reference_name(aligned_read.reference_id)
					read2_chr = samfile.get_reference_name(aligned_read.next_reference_id)

					if read1_chr.startswith(chr_prefix):
						if read2_chr.startswith(chr_prefix):
							aligned_read.reference_id = lookup_refid(ex_header, read1_chr)
							EX.write(aligned_read)
							fly_reads += 1
						else:
							BO.write(aligned_read)
							diff_genome += 1
					else:
						if read2_chr.startswith(chr_prefix):
							#aliged_read.reference_id = lookup_refid(ex_header, read2_chr)
							BO.write(aligned_read)
							diff_genome += 1
						else:
							aligned_read.reference_id = lookup_refid(human_header, read1_chr)
							HU.write(aligned_read)
							human_reads += 1
				else:
					read_chr = samfile.get_reference_name(aligned_read.reference_id)
					if read_chr.startswith(chr_prefix):
						aligned_read.reference_id = lookup_refid(ex_header, read_chr)
						EX.write(aligned_read)
						fly_reads += 1
					else:
						aligned_read.reference_id = lookup_refid(human_header, read_chr)
						HU.write(aligned_read)
						human_reads += 1
	except StopIteration:
		logging.info("Done")

	HU.close()
	EX.close()
	BO.close()
	NE.close()

	for i in ((outfile + '_human.bam'), (outfile + '_exogenous.bam'), (outfile + '_both.bam'), (outfile + '_neither.bam')):
		sort_bam(bamfile = i, threads = threads)


	print ("Sample\tn_unmapped\tn_qcFail\tn_duplicate\tn_secondary\tn_low_mapq\tn_both\tn_sample\tn_exogenous", file = OUT)
	print ("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (os.path.basename(bam_file),unmapped_reads, qcfail_reads, duplicate_reads, secondary_reads, low_maq, diff_genome, human_reads, fly_reads), file=OUT)
	OUT.close()

def main():

	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)

	parser.add_option("-i",action="store",type="string",dest="bam_file",help="BAM file of the composite genome (such as human + fly)")
	parser.add_option("-o","--output", action="store",type="string",dest="out_prefix",help="Output prefix. The original BAM file will be split into four BAM files: 'prefix_human.bam', 'prefix_exogenous.bam', 'prefix_both.bam', 'prefix_neither.bam'. ")
	parser.add_option("-p","--exo-prefix",action="store",type="string",dest="chr_prefix",default = 'dm6_', help="Prefix added to the exogenous chromosome IDs. For example. 'chr2L' -> 'dm6_chr2L'. default=%default")
	parser.add_option("-q", "--mapq",action="store",type="int",dest="map_qual",default=30, help="Mapping quality (phred scaled) threshold. Alignments with mapping quality score lower than this will be assigned to 'prefix_neither.bam'. default=%default")
	parser.add_option("--threads",action="store",type="int",dest="n_thread",default=1, help="Number of threads to use for BAM sorting. default=%default")
	(options,args)=parser.parse_args()

	if not (options.bam_file and options.out_prefix):
		parser.print_help()
		sys.exit()
	divided_bam(bam_file = options.bam_file, outfile = options.out_prefix, q_cut= options.map_qual, chr_prefix= options.chr_prefix, threads = options.n_thread)

if __name__=='__main__':
	main()

#sort and index bam files
