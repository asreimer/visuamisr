# -*- coding: utf-8 -*-
# Copyright (C) 2019  Ashton S. Reimer
# Full license can be found in LICENSE.txt
import os
import subprocess
from setuptools import find_packages, setup

requirements = ['pathlib2','cartopy','numpy','matplotlib','h5py','cython']

try:
    subprocess.call(['conda','install',' '.join(requirements)])
    requirements = []
except Exception:
    pass


readme = os.path.join(os.path.dirname(__file__), 'README.md')
with open(readme,'r') as f:
    readme_txt = f.readlines()
readme_txt = '\n'.join(readme_txt)
    

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
      install_requires=requirements,
      setup_requires=requirements,
      version="1.1",
      description="A library of data plotting utilities for visualizing processed Advanced Modular Incoherent Scatter Radar (AMISR) data.",
      author="Ashton S. Reimer",
      author_email="ashtonsethreimer@gmail.com",
      url="https://github.com/asreimer/pyAMISR",
      download_url="https://github.com/asreimer/pyAMISR",
      packages=find_packages(),
      long_description=readme_txt,
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
