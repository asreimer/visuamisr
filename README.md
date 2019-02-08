visuamisr
=========
A library of data plotting utilities for visualizing processed Advance Modular Incoherent Scatter Radar (AMISR) data.

Tested and working on both Python 2.7 and Python 3.7.

Install
=======
First clone this repository::

    git clone https://github.com/asreimer/visuamisr.git

Next run the `setup.py` file::

    python setup.py develop

or with pip::

    pip install .


Dependencies
------------
Makes use of the python wrapper for aacgm and the `mapObj` object in `davitpy`.


Usage
=====

First, you will need some data! I grabbed a fitted file from https://amisr.com/database, specifically, some PFISR Themis36 data from 2 March, 2016: https://amisr.com/database/61/experiment/20160302.001/3/2. The 20160302.001_lp_1min-fitcal.h5 file is 192 MB in size.

Beam Plot in Polar Coordinates
------------------------------
A visualization of the beam pattern used by the radar can be made in polar coordinates::

    import visuamisr
    isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
    isr.plot_polar_beam_pattern(min_elevation=10)

RTI Plotting
------------
Range Time Intensity (RTI) plots are a great way to visualize the data products of an incoherent scatter radar.
To make an RTI plot in `visuamisr` for one beam of data::

    import visuamisr
    from datetime import datetime
    isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
    isr.rti(['density','Te','Ti','velocity'],
            time_lim=[datetime(2016,3,2,6,0),datetime(2016,3,2,17)],
            ylim=[100,500],bmnum=10)

Profile Plotting
----------------
The altitude profile of various parameters can be plotted. For example::

    import visuamisr
    from datetime import datetime
    isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
    isr.profile_plot(['density','Te','Ti','velocity'],
                     datetime(2016,3,2,14,55),bmnum=10,
                     param_lim=[[10**10,10**12],[0,5000],[0,4000],
                                [-1000,1000]],use_range=True)

3D Beam Plotting
----------------
A 3 dimensional plot of the beams of the radar colour coded by a plasma parameter can be made::

    import visuamisr
    from datetime import datetime
    isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
    isr.plot_beams3d('density',datetime(2016,3,2,14,55),sym_size=5,clim=[10,12])
