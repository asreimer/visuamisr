#pyAMISR
A library of data plotting utilities for visualizing processed Advance Modular Incoherent Scatter Radar (AMISR) data. 

#Installation Instructions
First clone this repository using `git`:

    git clone https://github.com/asreimer/pyAMISR.git
Next run the `setup.py` file. If you would like to install `pyAMISR` to the system python site-packages directory, then run

    python setup.py install
as root, or run

    python setup.py install --user
as a regular user to install `pyAMISR` to the user python site-packages directory.

###Dependencies
This code depends on `numpy`, `matplotlib`, and `h5py`.


#Usage
###Beam Plot in Polar Coordinates
A visualization of the beam pattern used by the radar can be made in polar coordinates:

    import pyAMISR
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.plotPolarBeamPattern(maxZen=85)

###RTI Plotting
Range Time Intensity plots are a great way to visualize the data products of an incoherent scatter radar. To make an RTI plot in `pyAMISR`:

    import pyAMISR
    from datetime import datetime
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.rti(['density','Te','Ti','velocity'],
            timeLim=[datetime(2012,11,24,6,0),datetime(2012,11,24,7)],
            yLim=[100,500],bmnum=33)

###Profile Plotting
The altitude profile of various parameters can be plotted. For example:

    import pyAMISR
    from datetime import datetime
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.profilePlot(['density','Te','Ti','velocity'],
                    datetime(2012,11,24,6,5,0),bmnum=40,
                    paramLim=[[10**10,10**12],[0,5000],[0,4000],
                              [-1000,1000]],rang=True)

###3D Beam Plotting
A 3 dimensional plot of the beams of the radar colour coded by a plasma parameter can be made:

    import pyAMISR
    from datetime import datetime
    isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
    isr.plotBeams3D('density',datetime(2012,11,24,6,40),symSize=5, cLim=[10,12])
