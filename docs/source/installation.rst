
Prerequisites
--------------

These programs must be installed and callable from the command line:

- `Python 3 <https://www.python.org/downloads/>`_
- `Samtools <http://www.htslib.org/>`_
- `bowtie2 <https://github.com/BenLangmead/bowtie2>`_ (Only needed if your input is fastq file)


Python Dependencies
--------------------

- `MACS2 <https://pypi.org/project/MACS2/>`_
- `pysam <https://pypi.org/project/pysam/>`_
- `pyBigWig <https://github.com/deeptools/pyBigWig>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <https://www.scipy.org/>`_

.. note::
   These packages will be installed automatically if you use `pip3 <https://pip.pypa.io/en/stable/installing/>`_.

Install Spiker
---------------
Use pip3 to install spiker from `PyPI <https://pypi.org/project/spiker/>`_ or `github <https://github.com/liguowang/spiker>`_. ::

 $ pip3 install spiker
 or 
 $ pip3 install git+https://github.com/liguowang/spiker.git

Upgrade Spiker
---------------
::

 $ pip3 install spiker --upgrade	
 
Uninstall 
----------
::

 $ pip3 uninstall spiker


