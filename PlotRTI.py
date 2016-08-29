#!/usr/bin/env python
"""
demo RTI plot

./PlotRTI.py  myfile.h5
"""
import pyAMISR
from datetime import datetime
import sys

isr = pyAMISR.analyze(sys.argv[1])

isr.rti(['density','Te','Ti','velocity'],
        timeLim=None, #[datetime(2012,11,24,6,0),datetime(2012,11,24,7)],
        yLim=[100,500],
        #bmnum=33,
        )