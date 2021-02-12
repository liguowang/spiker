import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for spiker  -- Analysis workflow for Spike-in ChIP-seq data
"""

def main():
	setup(  name = "spiker",
			version = "1.0.0",
			py_modules = [ 'psyco_full' ],
			python_requires='>=3.5',
			packages = find_packages( 'lib' ),
			package_dir = { '': 'lib' },
			package_data = { '': ['*.ps'] },
			scripts = glob.glob( "bin/*.py"),
			ext_modules = [],
			test_suite = 'nose.collector',
			setup_requires = ['nose>=0.10.4'],
			author = "Liguo Wang",
			author_email ="wangliguo78@gmail.com",
			platforms = ['Linux','MacOS'],
			requires = ['cython (>=0.17)'],
			install_requires = ['numpy','scipy', 'pysam', 'deeptools', 'macs2',], 
			description = "spiker (Analysis workflow for ChIP-seq data with spike-in)",
			long_description = "Spiker is a tool to ananzlye data generated from ChIP-seq experiments with exogenous spike-in chromatin as internal control",
			license='MIT License',
			url = "",
			zip_safe = False,
			dependency_links = [],
			classifiers=[
				'Development Status :: 5 - Production/Stable',
				'Environment :: Console',
				'Intended Audience :: Science/Research',
				'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: POSIX',
				'Programming Language :: Python',
				'Topic :: Scientific/Engineering :: Bio-Informatics',
			],
			
			keywords='ChIP-seq, ChIPseq, spike-in, spikein, chromatin, peak calling',
             )


if __name__ == "__main__":
	main()
