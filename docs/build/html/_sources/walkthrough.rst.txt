.. role:: raw-math(raw)
    :format: latex html

Prepare FASTQ files
===================
Usually, the FASTQ files are delivered to you from the sequencing core. In this tutorial, we will use the H3K79me2 ChIP-seq data published by `Orlando et al (Cell Reports, 2014) <https://pubmed.ncbi.nlm.nih.gov/25437568/>`_. In this study, Drosophila S2 cells were added to human Jurkat cells at the ratio of 1:2.


+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+
| Sample (ChIP)      | SRR_accession                                                               | Sample (WCE)       | SRR_accession                                                                |
+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+
| Jurkat_K79_0%_R1   | `SRR1536557   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536557>`_ | Jurkat_WCE_0%_R1   | `SRR1584489    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584489>`_ |
+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+
| Jurkat_K79_25%_R1  | `SRR1536558   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536558>`_ | Jurkat_WCE_25%_R1  | `SRR1584490    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584490>`_ |
+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+
| Jurkat_K79_50%_R1  | `SRR1536559   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536559>`_ | Jurkat_WCE_50%_R1  | `SRR1584491    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584491>`_ |
+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+
| Jurkat_K79_75%_R1  | `SRR1536560   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536560>`_ | Jurkat_WCE_75%_R1  | `SRR1584492    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584492>`_ |
+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+
| Jurkat_K79_100%_R1 | `SRR1536561   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536561>`_ | Jurkat_WCE_100%_R1 | `SRR1584493    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584493>`_ |
+--------------------+-----------------------------------------------------------------------------+--------------------+------------------------------------------------------------------------------+

After download the SRA files to your hard drive, you will need the :code:`fastq-dump` command from the `SRA Tookit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc>`_ to convert the SRA archive files into FASTQ format. We used this command to convert single-end sequencing data:
::

 $ fastq-dump SRA_file1 SRA_file2 SRA_file3

.. note::
   Alternatively, you can download FASTQ files **directly** from the `ENA (Enropean Nucleotide Archive) <https://www.ebi.ac.uk/ena/browser/home>`_ by searching SRR accession numbers listed in the above table. 


Build bowtie2 index files
=========================
To calculate how many reads in the FASTQ files are drived from the Drosophila S2 cells, we could map all reads to the **composite reference genome** (i.e., human + Drosophila). In this tutorial, we will use `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_, other short reads aligners such as `BWA <http://bio-bwa.sourceforge.net/>`_ also work fine. To use `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ to map reads to the 
**composite reference genome**, we need to create bowtie2 index files first. 


Download the human and Drosophila reference genome sequences::

 $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
 $ wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
 $ gunzip hg38.fa.gz
 $ gunzip dm6.fa.gz

Rename chromosome IDs of the Drosophila genome::

 $ sed 's/chr/dm6_chr/g' dm6.fa  > dm6_renameID.fa

Concatenate the human and Drosophila genome files::

 $ cat hg38.fa  dm6_renameID.fa > hg38_dm6.fa

Run :code:`bowtie2-build` to build index files. "hg38_dm6" is the prefix of index files::

 $ bowtie2-build  hg38_dm6.fa  hg38_dm6

Clean up::

 $ rm hg38.fa dm6.fa dm6_renameID.fa

.. note::
   We have :ref:`bowtie2_index` files for **hg38 + dm6**, **hg19 + dm6**, **mm9 + dm6**, **mm10 + dm6**, **mm39 + dm6** available for download.

   
Map reads to the composite reference genome 
===========================================

Use `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ to map all reads to the **composite reference genome** (i.e., hg38_dm6.fa). Sort and index the BAM files using `Samtools <http://www.htslib.org/>`_. Below is an example to process one sample::

 bowtie2 --threads 4 -x ./database/GRCh38_dm6  -U Jurkat_K79_00p_Rep1_SRR1536557.fastq.gz  | samtools view -Sbh - > Jurkat_K79_00p_Rep1_SRR1536557.bam
 samtools sort -@ 4 Jurkat_K79_00p_Rep1_SRR1536557.bam > Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam
 samtools index -@ 4 Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam

.. note::
   We have pre-aligned BAM files available from :ref:`dataset` for download. 

Calculate scaling factors
=========================

To remove experimental bias and accurately quantify ChIP signal changes, the amount of exogenous Drosophila chromatin in different samples needs to be normalized to the same level. After sequencing, this is equivalent to normalize Drosophila reads to the same amount. 

We provide :code:`split_bam.py` to split the BAM file into 4 sub-BAM files. 

 * _both.bam : contains reads mapped to both human and Drosophila genomes
 * _exogenous.bam : contains reads only map to the Drosophila genomes
 * _human.bam : contains reads only map to the human genomes
 * _neither.bam : contains qcfailed, low quality, unmapped, dupilcate reads etc. 
 
A report file describing the number of alignments in each sub-BAM is also generated. 

For example::
 
 $ split_bam.py --threads 8  -i Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam  -o Jurkat_K79_00p_Rep1
 Read the BAM file: Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam
 Done
 [bam_sort_core] merging from 0 files and 8 in-memory blocks...
 [bam_sort_core] merging from 0 files and 8 in-memory blocks...
 [bam_sort_core] merging from 0 files and 8 in-memory blocks..

 $ cat Jurkat_K79_00p_Rep1.report.txt

 Sample	n_unmapped	n_qcFail	n_duplicate	n_secondary	n_low_mapq	n_both	n_sample	n_exogenous
 Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam	719553	0	0	0	5355549	0	20741684	5207395


Mapping statistics and the derived **scaling factors (SF)** are summerized in the table below (the number of Drosophila reads in each sample was normalized to 1 million).

+---------------------------------+----------------+----------------+-------------+------------+------------+
| Sample                          | Unmapped reads | Low_mapq_reads | Human_reads | Fly_reads  | SF         |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_K79_00p_Rep1_SRR1536557  | 719,553        | 5,355,549      | 20,741,684  | 5,207,395  | 0.1920346  |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_K79_25p_Rep1_SRR1536558  | 2,177,315      | 5,367,333      | 19,282,608  | 6,035,961  | 0.1656737  |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_K79_50p_Rep1_SRR1536559  | 2,919,492      | 5,233,017      | 17,235,944  | 6,931,602  | 0.14426679 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_K79_75p_Rep1_SRR1536560  | 1,627,494      | 4,300,454      | 11,731,507  | 8,291,870  | 0.12060006 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_K79_100p_Rep1_SRR1536561 | 4,526,721      | 5,109,981      | 7,893,377   | 14,372,323 | 0.06957817 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_WCE_00p_Rep1.SRR1584489  | 2,344,210      | 3,914,206      | 11,126,768  | 844,395    | 1.18427987 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_WCE_25p_Rep1.SRR1584490  | 9,982,131      | 1,056,402      | 2,771,236   | 164,052    | 6.09562822 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_WCE_50p_Rep1.SRR1584491  | 4,398,297      | 1,750,494      | 4,431,034   | 267,211    | 3.74236091 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_WCE_75p_Rep1.SRR1584492  | 7,317,599      | 1,620,854      | 4,029,096   | 179,641    | 5.56665795 |
+---------------------------------+----------------+----------------+-------------+------------+------------+
| Jurkat_WCE_100p_Rep1.SRR1584493 | 7,609,117      | 2,336,349      | 5,197,374   | 344,665    | 2.901368   |
+---------------------------------+----------------+----------------+-------------+------------+------------+

.. Note::
   You can also use K-mer based, alignment-free methods to calculate scaling factor. For example, you can use `xenome <https://github.com/data61/gossamer/blob/master/docs/xenome.md>`_ to separate human and fly reads from your FASTQ files.

Run Spiker
==========

**Spiker** is specifically designed to analyze ChIP-seq data with spike-in control. The input files can be single-end or
paired-end `FASTQ <https://en.wikipedia.org/wiki/FASTQ_format#>`_ files or pre-aligned `BAM <https://genome.ucsc.edu/goldenPath/help/bam.html>`_ files. It supports both narrow peak (such as H3K27ac) and broad peak analysis (such as H3K79me2 in this tutorial). 

In below, we used "Jurkat_K79_00p_Rep1_SRR1536557_human.sorted.bam" (ChIP) and "Jurkat_WCE_00p_Rep1.SRR1584489_human.sorted.bam" (control) as inputs to Spiker::


 $ spiker.py --broad --spikeIn --csf 1.18 --tsf 0.19 --bw -t Jurkat_K79_00p_Rep1_SRR1536557_human.sorted.bam -c Jurkat_WCE_00p_Rep1.SRR1584489_human.sorted.bam  -o Jurkat_K79_00p_Rep1_SRR1536557

 2021-02-18 12:56:20 [INFO]  Running ChIP-seq workflow ...
 2021-02-18 12:56:20 [INFO]  Input BAM files. Skip Step_1 (reads mapping).
 2021-02-18 12:56:20 [INFO]  Step_2: Extract information from BAM files ...
 2021-02-18 12:56:20 [INFO]  Step_2.1: Check sequencing layout from BAM file ...
 2021-02-18 12:56:22 [INFO]  	The layout of Jurkat_K79_00p_Rep1_SRR1536557_human.sorted.bam is SE
 2021-02-18 12:56:23 [INFO]  	The layout of Jurkat_WCE_00p_Rep1.SRR1584489_human.sorted.bam is SE
 2021-02-18 12:56:23 [INFO]  Step_2.2: Extract information from ChIP BAM files ...
 2021-02-18 12:56:25 [INFO]  	The total mapped reads in ChIP sample is: 20741684
 2021-02-18 12:56:25 [INFO]  Step_2.3: Extract information from control BAM files ...
 2021-02-18 12:56:26 [INFO]  	The total mapped reads in control sample is: 11126768
 2021-02-18 12:56:26 [INFO]  Step_3: Deduplication ...
 2021-02-18 12:56:26 [INFO]  Step_3.1: Remove duplidate reads from ChIP BAM file ...
 ...

Output files

.. glossary::

	*.gappedPeak
	  peak file (in `gappedPeak <https://genome.ucsc.edu/FAQ/FAQformat.html#format14>`_) called by Spiker. 

	*.control.pileup.max.bdg
	  `bedGraph <https://genome.ucsc.edu/goldenPath/help/bedgraph.html>`_ file of **control pileup**

	*.control.pileup.max.bigWig
	  `bigWig <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`_ file of **control pileup** 

	*.treat.pileup.bdg
	  `bedGraph <https://genome.ucsc.edu/goldenPath/help/bedgraph.html>`_ file of **ChIP sample pileup** (raw)

	*.treat.pileup.SpikeIn_scaled.bdg
	  `bedGraph <https://genome.ucsc.edu/goldenPath/help/bedgraph.html>`_ file of **ChIP sample pileup** (Scaled to spike-in control)

	*.treat.pileup.SpikeIn_scaled.bigWig
	  `bigWig <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`_ file of **ChIP sample pileup** (Scaled to spike-in control)

