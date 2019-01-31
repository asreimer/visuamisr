#!/usr/bin/env python
""" visamisr is a data visualization tool for AMISR data
It provides:
- A data reading utilty
- Range Time Intensity plotting
- Profile plotting
- 3D beam plotting
The full license can be found in LICENSE.txt
"""

import os
import subprocess
from setuptools import find_packages, setup

REQUIREMENTS = ['pathlib2', 'numpy', 'matplotlib', 'h5py']

try:
    subprocess.call(['conda', 'install', ' '.join(REQUIREMENTS)])
    REQUIREMENTS = []
except subprocess.builtins.FileNotFoundError:
    pass


README = os.path.join(os.path.dirname(__file__), 'README.md')
with open(README, 'r') as f:
    READMETXT = f.readlines()
READMETXT = '\n'.join(READMETXT)


DESC = "A library of data plotting utilities for visualizing processed "
DESC += "Advanced Modular Incoherent Scatter Radar (AMISR) data."

#############################################################################
# First, check to make sure we are executing
# 'python setup.py install' from the same directory
# as setup.py (root directory)
#############################################################################
PATH = os.getcwd()
assert('setup.py' in os.listdir(PATH)), \
       "You must execute 'python setup.py install' from within the \
repo root directory."


#############################################################################
# Now execute the setup
#############################################################################
setup(name='visamisr',
      install_requires=REQUIREMENTS,
      setup_requires=REQUIREMENTS,
      version="2.0.0",
      description=DESC,
      author="Ashton S. Reimer",
      author_email="ashtonsethreimer@gmail.com",
      url="https://github.com/asreimer/visamisr",
      download_url="https://github.com/asreimer/visamisr",
      packages=find_packages(),
      long_description=READMETXT,
      zip_safe=False,
      py_modules=['visamisr'],
      classifiers=["Development Status :: 2.0.0 - Release",
                   "Topic :: Scientific/Engineering",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Natural Language :: English",
                   "Programming Language :: Python",
                  ],
      )
