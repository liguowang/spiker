.. _bowtie2_index:

Bowtie2 index files
-------------------

We first download the Reference genome sequences for Human, Mouse, and Drosophila from `UCSC <http://hgdownload.soe.ucsc.edu/downloads.html>`_. We then build the `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ index files for **human + Drosophila** and **mouse + Drosophila** composite genomes (listed in the table below). 

Index files for Human + Drosophila
------------------------------------------

+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Composite genome                       | Download link                                                                                             | Size   | md5sum                           |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Human (GRCh37/hg19) + Drosophila (dm6) | `GRCh37_dm6.tar.gz <https://de.cyverse.org/dl/d/923CED1F-23C8-4637-915B-53D1BA571207/GRCh37_dm6.tar.gz>`_ | 4.9 GB | 962d0356703462f610cde06e46bde9e4 |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Human (GRCh38/hg38) + Drosophila (dm6) | `GRCh38_dm6.tar.gz <https://de.cyverse.org/dl/d/6406D535-893A-4509-AD7F-4FE6D604ADE3/GRCh38_dm6.tar.gz>`_ | 5.1 GB | efb62dac65aafae0ca8d300525cae7c6 |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+


Index files for Mouse + Drosophila
------------------------------------------

+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Composite genome                       | Download link                                                                                             | Size   | md5sum                           |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Mouse (NCBI37/mm9) + Drosophila (dm6)  | `mm9_dm6.tar.gz <https://de.cyverse.org/dl/d/AD2188D4-E5D6-43AD-A10C-C39C30383B27/mm9_dm6.tar.gz>`_       | 4.3 GB | 885d11585ed5d407f9f37eafd8fff106 |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Mouse (GRCm38/mm10) + Drosophila (dm6) | `mm10_dm6.tar.gz <https://de.cyverse.org/dl/d/008771A7-19D2-48D3-BD72-C2038CF22951/mm10_dm6.tar.gz>`_     | 4.5 GB | 751753ff50bb2c228e213a8873e92943 |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+
| Mouse (GRCm39/mm39) + Drosophila (dm6) | `mm39_dm6.tar.gz <https://de.cyverse.org/dl/d/654495F8-FD14-4C8C-9F54-49B0B6E60C6A/mm39_dm6.tar.gz>`_     | 4.5 GB | e1cafa46619ad1a1bfd98e00a635e853 |
+----------------------------------------+-----------------------------------------------------------------------------------------------------------+--------+----------------------------------+

After download, use :code:`md5sum *.tar.gz` (Linux) or :code:`md5 *.tar.gz` (Mac OS X) to verify the file integrity and then use :code:`tar -zxvf *.tar.gz` to uncompress.

.. Note::
   Because the Drosophila genome also has the "chrX", "chrY" and, "chrM". In the composite genome, the chromosome IDs of Drosophila were modified by adding the prefix "dm6_". 
