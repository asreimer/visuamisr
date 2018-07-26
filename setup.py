# -*- coding: utf-8 -*-
# Copyright (C) 2013  Ashton S. Reimer
# Full license can be found in LICENSE.txt
import os
from setuptools import find_packages, setup
import subprocess

requirements = ['pathlib2','cartopy','numpy','matplotlib','h5py']

try:
    subprocess.call(['conda','install',' '.join(requirements)])
    requirements = []
except Exception:
    pass


#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root directory)
#############################################################################
path = os.getcwd()
assert('setup.py' in os.listdir(path)), \
       "You must execute 'python setup.py install' from within the \
repo root directory."


#############################################################################
# Now execute the setup
#############################################################################
setup(name='pyAMISR',
      install_requires=requirements
      version="1.1",
      description="A library of data plotting utilities for visualizing processed Advanced Modular Incoherent Scatter Radar (AMISR) data.",
      author="Ashton S. Reimer",
      author_email="ashtonsethreimer@gmail.com",
      url="https://github.com/asreimer/pyAMISR",
      download_url="https://github.com/asreimer/pyAMISR",
      packages=find_packages(),
      long_description=read('README.rst'),
      zip_safe=False,
      py_modules=['pyAMISR'],
      classifiers=[
            "Development Status :: 1.1 - Release",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Programming Language :: Python"
            ],
      )
