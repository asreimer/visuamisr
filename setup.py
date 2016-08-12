# -*- coding: utf-8 -*-
# Copyright (C) 2013  Ashton S. Reimer
# Full license can be found in LICENSE.txt
import os
from setuptools import find_packages, setup
import subprocess

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


try:
    subprocess.call(['conda','install','--file','requirements.txt'])
    req = []
except Exception:
    req = read('requirements.txt')


#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root directory)
#############################################################################
path = os.getcwd()
assert('setup.py' in os.listdir(path)), \
       "You must execute 'python setup.py install' from within the \
davitpy root directory."


#############################################################################
# Now execute the setup
#############################################################################
setup(name='pyAMISR',
      install_requires=req,
      version="1.0",
      description="A library of data plotting utilities for visualizing processed Advanced Modular Incoherent Scatter Radar (AMISR) data.",
      author="VT SuperDARN Lab and friends",
      author_email="ashtonsethreimer@gmail.com",
      url="https://github.com/asreimer/pyAMISR",
      download_url="https://github.com/asreimer/pyAMISR",
      packages=find_packages(),
      long_description=read('README.rst'),
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
