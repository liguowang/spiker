
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:11:27 2019

@author: m102324
"""
import sys,os
import logging
import shutil
import subprocess


def bowtie2_map(bt2_prefix, fq1, fq2, outbam, nthread = 8):
	'''
	Run bowtie2, samtools sort, and samtools index commands

	parameters
	----------

	bt2_prefix : str
		Index filename prefix (minus trailing .X.bt2) of bowtie2
		eg: /research/bsi/development/m102324/database/TCGA_database/bowtie_index/GRCh38

	fq1 : str
		fastq file of read1

	fq2 : str
		fastq file of read2. No need for single-end sequenicng data.

	outbam : str
		Output BAM file name

	nthread : int
		 Number of threads

	print_cmd : bool
		if set to True, return the "bowtie2" command line and exit. If set to False
		run the "bowtie2" command line.
	'''

	logging.info("\tChecking input files ...")
	# find bowtie2 command
	bowtie2_cmd = shutil.which("bowtie2")
	if bowtie2_cmd is None:
		logging.error("\tCannot find the \"bowtie2\" command!")
		sys.exit()
	logging.debug("\tFound bowtie2 command:\"%s\"" %  bowtie2_cmd)

	# find samtools command
	samtools_cmd = shutil.which("samtools")
	if samtools_cmd is None:
		logging.error("\tCannot find the \"samtools\" command!")
		sys.exit()
	logging.debug("\tFound samtools command:\"%s\"" %  samtools_cmd)



	if outbam.endswith('.bam'):
		base_bam = outbam
		sorted_bam = outbam.replace('.bam','.sorted.bam')
		sorted_bam_bai = outbam.replace('.bam','.sorted.bam.bai')
	else:
		base_bam = outbam + '.bam'
		sorted_bam = outbam + '.sorted.bam'
		sorted_bam_bai = outbam + '.sorted.bam.bai'

	if os.path.exists(base_bam) and os.path.getsize(base_bam) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip bowtie2 mapping." % base_bam)
		pass
	elif os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip bowtie2 mapping." % sorted_bam)
		pass
	else:
		fq_files = []

		for ff in [(bt2_prefix  + i) for i in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]]:
			if not os.path.exists(ff):
				logging.error("\tBowtie2 index file: \"%s\" does not exist!" % ff)
				sys.exit()

		if not os.path.exists(fq1):
			logging.error("\tCannot find \"%s\"" % fq1)
			sys.exit()
		else:
			fq_files.append(os.path.abspath(fq1))

		if (fq2 is None) or (not os.path.exists(fq2)):
			logging.warning("\tCannot find \"%s\", use single-end read mapping mode." % fq2)
			pass


		# run bowtie2
		if len(fq_files) == 2:
			cmd = ' '.join([bowtie2_cmd, '--threads', str(nthread), '-x', bt2_prefix, '-1', fq_files[0], '-2', fq_files[1], '| samtools view -Sbh -o', base_bam, '-'])
		elif len(fq_files) == 1:
			cmd = ' '.join([bowtie2_cmd, '--threads', str(nthread), '-x', bt2_prefix, '-U', fq_files[0], '| samtools view -Sbh -o', base_bam, '-'])
		else:
			logging.error("\tNo fastq files were provided.")
			sys.exit()

		logging.info("\tRuning bowtie2 ...")
		logging.info("\tCommand: %s" % cmd)
		subprocess.call(cmd, shell=True)

	# sort BAM file
	if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip BAM sorting." % sorted_bam)
		pass
	else:
		samtools_sort = "%s sort -@ %d %s > %s" % (samtools_cmd, nthread, base_bam, sorted_bam)
		logging.info("\tSorting BAM ...")
		logging.info("\tCommand: %s" %  samtools_sort)
		subprocess.call(samtools_sort, shell=True)

	# index BAM file
	if os.path.exists(sorted_bam_bai) and os.path.getsize(sorted_bam_bai) > 0:
		logging.warning("\t\"%s\" exists and non-empty, skip BAM indexing." % sorted_bam_bai)
		pass
	else:
		samtools_index = "%s index %s" % (samtools_cmd, sorted_bam)
		logging.info("\tIndexing BAM ...")
		logging.info("\tCommand: %s" %  samtools_index)
		subprocess.call(samtools_index, shell=True)


if __name__=='__main__':
	if len(sys.argv) < 5:
		print ('bowtie2_map.py  bt2_prefix  fq1  fq2  outbam  nthread', file = sys.stderr)
		sys.exit(0)
	bowtie2_map(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])