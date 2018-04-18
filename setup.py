#!/usr/bin/env python
from setuptools import setup
import sysconfig
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
from glob import glob
import os
version = '2.1.1'

# Since Conda started using their own compilers the include dir may not be found when trying to install within an
# active conda environment.
# Here we explicitly set the include path for the current environment.
extra_include_path = os.path.join(sysconfig.get_path("data"),"include")
extra_compile_args = "-I{}".format(extra_include_path)


module1 = Extension("sampletau",
                             sources=["sampletau/sampletau.pyx", "sampletau/c_sample_tau.c"],
                             extra_compile_args=['-O3', '-march=native', extra_compile_args],
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
      scripts=["bin/desman","desman/Variant_Filter.py", "desman/GeneAssign.py"]+glob("scripts/*"),
      include_package_data=True,
      zip_safe=False,
      cmdclass={'build_ext': build_ext},
      ext_modules=[module1],
      extra_compile_args=extra_compile_args,
      install_requires=['cython>=0.19.1',
                        'numpy>=1.7.1',
                        'scipy>=0.12.0',
                        'pandas>=0.11.0'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
