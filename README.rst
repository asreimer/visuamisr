=======
pyAMISR
=======
A library of data plotting utilities for visualizing processed Advance Modular Incoherent Scatter Radar (AMISR) data.

Install
=======
First clone this repository::

    git clone https://github.com/asreimer/pyAMISR.git

Next run the `setup.py` file::

    python setup.py develop


Dependencies
------------
Makes use of the python wrapper for aacgm and the `mapObj` object in `davitpy`.


Usage
=====

Beam Plot in Polar Coordinates
------------------------------
A visualization of the beam pattern used by the radar can be made in polar coordinates::

    import pyAMISR
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.plotPolarBeamPattern(maxZen=85)

RTI Plotting
------------
Range Time Intensity plots are a great way to visualize the data products of an incoherent scatter radar.
To make an RTI plot in `pyAMISR`::

    import pyAMISR
    from datetime import datetime
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.rti(['density','Te','Ti','velocity'],
            timeLim=[datetime(2012,11,24,6,0),datetime(2012,11,24,7)],
            yLim=[100,500],bmnum=33)

Profile Plotting
----------------
The altitude profile of various parameters can be plotted. For example::

    import pyAMISR
    from datetime import datetime
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.profilePlot(['density','Te','Ti','velocity'],
                    datetime(2012,11,24,6,5,0),bmnum=40,
                    paramLim=[[10**10,10**12],[0,5000],[0,4000],
                              [-1000,1000]],rang=True)

3D Beam Plotting
----------------
A 3 dimensional plot of the beams of the radar colour coded by a plasma parameter can be made::

    import pyAMISR
    from datetime import datetime
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.plotBeams3D('density',datetime(2012,11,24,6,40),symSize=5, cLim=[10,12])
