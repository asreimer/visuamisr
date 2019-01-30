# -*- coding: utf-8 -*-
# Copyright (C) 2013  Ashton S. Reimer
# Full license can be found in LICENSE.txt
"""
pyamisr
------

Module for plotting AMISR data

Modules
-------------------------------------------
read_data     Read data from HDF5 file
analyze       Plotting utilities
---------- --------------------------------

"""
try:
    from pathlib import Path
    Path().expanduser()
except (ImportError, AttributeError):
    from pathlib2 import Path

from .analyze import analyze, read_data

__version__ = '2.0.0'
