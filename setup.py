# -*- coding: utf-8 -*-
# Copyright (C) 2013  Ashton S. Reimer
# Full license can be found in LICENSE.txt
import os
# Need to use the enhanced version of distutils packaged with
# numpy so that we can compile fortran extensions
from setuptools import find_packages, setup

#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root directory)
#############################################################################
path = os.getcwd()
assert('setup.py' in os.listdir(path)), \
       "You must execute 'python setup.py install' from within the \
davitpy root directory."

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

#############################################################################
# Now execute the setup
#############################################################################
setup(name='pyAMISR',
      install_requires=['numpy','matplotlib','h5py'],
      version="1.0",
      description="A library of data plotting utilities for visualizing processed Advance Modular Incoherent Scatter Radar (AMISR) data.",
      author="VT SuperDARN Lab and friends",
      author_email="ashtonsethreimer@gmail.com",
      url="",
      download_url="https://github.com/asreimer/pyAMISR",
      packages=find_packages(),
      long_description=read('README.md'),
      zip_safe=False,
      py_modules=['pyAMISR'],
      classifiers=[
            "Development Status :: 1.0 - Release",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Programming Language :: Python"
            ],
      )
