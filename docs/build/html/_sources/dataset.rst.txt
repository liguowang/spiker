.. _dataset:

Test dataset
------------

We used the H3K79me2 ChIP-seq data published by `Orlando et al <https://pubmed.ncbi.nlm.nih.gov/25437568/>`_ as the testing dataset. In this study, human Jurkat cells were treated with either DMSO or EPZ (Dot1L inhibitor), then collected individually and mixed according to the following compositions by cell number.

   #. 100% DMSO to 0% EPZ
   #. 75% DMSO to 25% EPZ.
   #. 50% DMSO to 50% EPZ.
   #. 25% DMSO to 75% EPZ.
   #. 0% DMSO to 100% EPZ.

They tested whether traditional ChIP-seq analysis methods would reveal the decrease in human per-cell H3K79me2 occupancy and, if not, whether the addition of the exogenous Drosophila chromatin would allow detection of H3K79me2 removal. All data has been deposited into GEO with accession `GSE60104 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60104>`_.


H3K79me2 ChIP-seq data
-------------------------

+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| Sample                | SRR_accession                                                               | BAM (hg38 only)                                                                                                               | bigWig files                                                                                                                 | BAM (hg38 + dm6)                                                                                                       |
+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                             |                                                                                                                               |                                                                                                                              |                                                                                                                        |
|    Jurkat_K79_0%_R1   | `SRR1536557   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536557>`_ | `bam   <https://de.cyverse.org/dl/d/600F3474-7E67-4673-B998-DBE963FBD76D/`Jurkat_K79_00p_Rep1_SRR1536557_human.sorted.bam>`_  | `bigWig   <https://de.cyverse.org/dl/d/940C4990-D88D-446E-A7EA-7626045D1237/`K79_00p_R1.treat.hg38.SpikeIn_scaled.bigWig>`_  | `bam   <https://de.cyverse.org/dl/d/A9C49FED-220C-4996-927A-1C4F160F5C30/Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam>`_  |
+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                             |                                                                                                                               |                                                                                                                              |                                                                                                                        |
|    Jurkat_K79_25%_R1  | `SRR1536558   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536558>`_ | `bam   <https://de.cyverse.org/dl/d/B1211D05-229A-4C3B-8A9F-3B5818CF3F58/`Jurkat_K79_25p_Rep1_SRR1536558_human.sorted.bam>`_  | `bigWig   <https://de.cyverse.org/dl/d/DCB02791-A9A4-4789-AF59-B2F651C14A27/`K79_25p_R1.treat.hg38.SpikeIn_scaled.bigWig>`_  | `bam   <https://de.cyverse.org/dl/d/141AFB03-67BE-416B-B00D-3B8C1553F521/Jurkat_K79_25p_Rep1_SRR1536558.sorted.bam>`_  |
+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                             |                                                                                                                               |                                                                                                                              |                                                                                                                        |
|    Jurkat_K79_50%_R1  | `SRR1536559   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536559>`_ | `bam   <https://de.cyverse.org/dl/d/77B7D500-8ECD-487D-8A10-F04F58CDA788/`Jurkat_K79_50p_Rep1_SRR1536559_human.sorted.bam>`_  | `bigWig   <https://de.cyverse.org/dl/d/B2545A6C-5ED3-4C6E-AB9D-925374DD8DCC/`K79_50p_R1.treat.hg38.SpikeIn_scaled.bigWig>`_  | `bam   <https://de.cyverse.org/dl/d/3041045F-7FD5-4FFA-A1BB-ED71DAE77D38/Jurkat_K79_50p_Rep1_SRR1536559.sorted.bam>`_  |
+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                             |                                                                                                                               |                                                                                                                              |                                                                                                                        |
|    Jurkat_K79_75%_R1  | `SRR1536560   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536560>`_ | `bam   <https://de.cyverse.org/dl/d/38734B6A-09B1-41E5-98CD-11B3CEE87442/`Jurkat_K79_75p_Rep1_SRR1536560_human.sorted.bam>`_  | `bigWig   <https://de.cyverse.org/dl/d/FE2A56A4-0BF1-4130-B6C9-6C9A08150499/`K79_75p_R1.treat.hg38.SpikeIn_scaled.bigWig>`_  | `bam   <https://de.cyverse.org/dl/d/42CC6AF3-E6FE-4CD1-823F-653E1ABBFF0C/Jurkat_K79_75p_Rep1_SRR1536560.sorted.bam>`_  |
+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                             |                                                                                                                               |                                                                                                                              |                                                                                                                        |
|    Jurkat_K79_100%_R1 | `SRR1536561   <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1536561>`_ | `bam   <https://de.cyverse.org/dl/d/270059B9-93C3-4D70-82CA-CDE4F30FE245/`Jurkat_K79_100p_Rep1_SRR1536561_human.sorted.bam>`_ | `bigWig   <https://de.cyverse.org/dl/d/315F8E60-D1BE-4E2F-A678-971A4F1620E9/`K79_100p_R1.treat.hg38.SpikeIn_scaled.bigWig>`_ | `bam   <https://de.cyverse.org/dl/d/AAEE8B55-D5E6-47C6-9DCE-4B33AF944835/Jurkat_K79_100p_Rep1_SRR1536561.sorted.bam>`_ |
+-----------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+


WCE (whole cell extract) data
-----------------------------

+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+
| Sample                | SRR_accession                                                                | BAM (hg38 only)                                                                                                              |
+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                              |                                                                                                                              |
|    Jurkat_WCE_0%_R1   | `SRR1584489    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584489>`_ | `bam   <https://de.cyverse.org/dl/d/757F5C35-4644-46A5-9CCD-8835F5EF16FC/Jurkat_WCE_00p_Rep1.SRR1584489_human.sorted.bam>`_  |
+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                              |                                                                                                                              |
|    Jurkat_WCE_25%_R1  | `SRR1584490    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584490>`_ | `bam   <https://de.cyverse.org/dl/d/F4835E45-F70E-4D41-97BC-96B3323F575D/Jurkat_WCE_25p_Rep1.SRR1584490_human.sorted.bam>`_  |
+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                              |                                                                                                                              |
|    Jurkat_WCE_50%_R1  | `SRR1584491    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584491>`_ | `bam   <https://de.cyverse.org/dl/d/DD765972-2588-4E06-B3A8-1701D8310263/Jurkat_WCE_50p_Rep1.SRR1584491_human.sorted.bam>`_  |
+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                              |                                                                                                                              |
|    Jurkat_WCE_75%_R1  | `SRR1584492    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584492>`_ | `bam   <https://de.cyverse.org/dl/d/C99A0893-CE57-4CDB-8A3D-57C9709FC9A2/Jurkat_WCE_75p_Rep1.SRR1584492_human.sorted.bam>`_  |
+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+
|                       |                                                                              |                                                                                                                              |
|    Jurkat_WCE_100%_R1 | `SRR1584493    <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1584493>`_ | `bam   <https://de.cyverse.org/dl/d/BFF0E876-7797-4961-B310-F0AB8DB6356E/Jurkat_WCE_100p_Rep1.SRR1584493_human.sorted.bam>`_ |
+-----------------------+------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------+

MD5sum
-------
::

 # K79 BAM files with reads only mapped to the human reference genome (GRCh38/hg38) 
 MD5 (Jurkat_K79_00p_Rep1_SRR1536557_human.sorted.bam) = a4b8ac56f6ce4e1f9c0da2279ec3ea36
 MD5 (Jurkat_K79_25p_Rep1_SRR1536558_human.sorted.bam) = 56be607cb3ab886aa2fb033b4a8d37ab
 MD5 (Jurkat_K79_50p_Rep1_SRR1536559_human.sorted.bam) = 653170d551ce395e40b044cfa36fbaef
 MD5 (Jurkat_K79_75p_Rep1_SRR1536560_human.sorted.bam) = 23495e1422d9546ba496e69e2becd6f7
 MD5 (Jurkat_K79_100p_Rep1_SRR1536561_human.sorted.bam) = 2ebc83a02982a2be8deaca7b84651e44
 
 # K79 bigWig files generated by Spiker
 MD5 (K79_00p_R1.treat.hg38.SpikeIn_scaled.bigWig) = f13d2c6e228c87cf1d8262352f7e983c
 MD5 (K79_25p_R1.treat.hg38.SpikeIn_scaled.bigWig) = 0789024f52c35a3888254239dad6766a
 MD5 (K79_50p_R1.treat.hg38.SpikeIn_scaled.bigWig) = 83321afec8522f21c20baee3e772d85c
 MD5 (K79_75p_R1.treat.hg38.SpikeIn_scaled.bigWig) = bb64aafd851acb581742b581fa06cc5b
 MD5 (K79_100p_R1.treat.hg38.SpikeIn_scaled.bigWig) = aa2a38350a0dc07978b76297ea913dc3
 
 # K79 BAM files with reads mapped to the composite genome (hg38 + dm6)
 MD5 (Jurkat_K79_00p_Rep1_SRR1536557.sorted.bam) = 729b5391769c763eb1c85e5b2ee66af9
 MD5 (Jurkat_K79_25p_Rep1_SRR1536558.sorted.bam) = faffab3d3f9161304826c7f5bf84ee9e
 MD5 (Jurkat_K79_50p_Rep1_SRR1536559.sorted.bam) = 0c72e603c844fa8a38ca9fe8c2496a89
 MD5 (Jurkat_K79_75p_Rep1_SRR1536560.sorted.bam) = 9cac7b03ec65804ad3cb25143480348a
 MD5 (Jurkat_K79_100p_Rep1_SRR1536561.sorted.bam) = 05107197f676e8b898bee99ed001a61a

 # WCE BAM files with reads only mapped to the human reference genome (GRCh38/hg38)
 MD5 (Jurkat_WCE_00p_Rep1.SRR1584489_human.sorted.bam) = 561430e24b5edab177bfc4e832c5d8c0
 MD5 (Jurkat_WCE_100p_Rep1.SRR1584493_human.sorted.bam) = 8a8546927383afeda6c0e20f82807024
 MD5 (Jurkat_WCE_25p_Rep1.SRR1584490_human.sorted.bam) = 4396f3a31c718029d620dd770825df1e
 MD5 (Jurkat_WCE_50p_Rep1.SRR1584491_human.sorted.bam) = 96c5ae3cd9ad129167d2143fcb66606b
 MD5 (Jurkat_WCE_75p_Rep1.SRR1584492_human.sorted.bam) = 09b7ec2ecec4ccb7cce8392f88a17125


