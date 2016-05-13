# -*- coding: utf-8 -*-
# Copyright (C) 2013  Ashton S. Reimer
# Full license can be found in LICENSE.txt
"""
pyAMISR
------

Module for plotting AMISR data

Modules
-------------------------------------------
read_data     Read data from HDF5 file
analyze       Plotting utilities
---------- --------------------------------

"""
try:
	from analyze import *
except Exception, e:
    print 'Problem importing from analyze.py: ' + str(e)