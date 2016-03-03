#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

version = '0.1dev'



module1 = Extension("sampletau",
                             sources=["sampletau/sampletau.pyx", "sampletau/c_sample_tau.c"],
                             extra_compile_args=['-O3', '-march=native'],
                             libraries =['gsl',  'gslcblas'],
                             include_dirs=[numpy.get_include(), '/opt/local/include/'])

setup(name='desman',
      version=version,
      description="De novo Extraction of Strains from MetAgeNomes",
      long_description="""DESMAN is a pipeline that facilities the
      completely de novo extraction of both strain haplotypes
      and accessory genomes de novo from metagenome data""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Metagenomics Python',
      author='Christopher Quince',
      author_email='c.quince@warwick.ac.uk',
      maintainer='Christopher Quince',
      maintainer_email='c.quince@warwick.ac.uk',
      url='https://github.com/chrisquince/DESMAN',
      license='FreeBSD',
      packages=['desman'],
      scripts=["bin/desman"],
      include_package_data=True,
      zip_safe=False,
      cmdclass = {'build_ext': build_ext},
      ext_modules=[module1],
      install_requires=['cython>=0.19.1',
                        'numpy>=1.7.1',
                        'scipy>=0.12.0',
                        'pandas>=0.11.0'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
