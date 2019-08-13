# -*- coding: utf-8 -*-
# Copyright (C) 2013  Ashton S. Reimer
# Full license can be found in LICENSE.txt
"""
.. module:: analyze
   :synopsis: Used for basic visualization of ISR data. accepts data in
              hdf5 format only. This module read data from the hdf5 in
              tp a dictionary with the following keys where (B =
              # beams, N = # time steps, R = # range gates):
              'az' - B length vector of azimuthal angles in degrees for
                     each ISR beam
              'el' - B length vector of elevation angles in degrees
                     for each ISR beam
              'site_latitude' - the ISR site geographic latitude
              'site_longitude' - the ISR site geographic longitude
              'site_altitude' - the ISR site altitude above sea level
                                in meters
              'times' - Nx2 array of datetime describing the start and
                        end times for an ISR measurement
              'ave_times' - N array of datetime describing the average
                           time for each ISR measurement
              'altitude' - BxR array of heights above sea level in
                           meters for each range cell on each beam
              'range' - BxR array of range along beam direction in
                        meters for each range cell on each beam
              'density' - NxBxR electron density (not log10)
              'edensity' - NxBxR error in electron density (not log10)
              'density_uncor' - NxBxR uncorrected electron density (not log10)
              'edensity_uncor' - NxBxR error in uncorrected electron density (not log10)
              'Te' - NxBxR electron temperature in Kelvin
              'eTe' - NxBxR error in electron density
              'Ti' - NxBxR ion temperature in Kelvin
              'eTi' - NxBxR error in ion temperature
              'vel' - NxBxR line of sight velocity in m/s
              'evel' - NxBxR error in line of sight velocity
              'latitude' - BxR geographic latitude of each range
                                  cell
              'longitude' - BxR geographic longitude of each
                                   range cell
              'babs' - BxR magnitude of magnetic field in each
                              range cell
              'kvec' - Bx3 k-vector of each beam in local North-
                              East-Up coordinate system.

"""

from __future__ import print_function

import calendar
from datetime import datetime

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cmx
from matplotlib import pyplot, colors, dates
from mpl_toolkits.mplot3d import Axes3D

from . import Path

# Dictionary of site code number/radar name pairs
SITECODES = {91: 'RISR-N', 92: 'RISR-C', 61: 'PFISR'}


# A function for reading fitted data from SRI format hdf5 files
def read_data(filepath):
    """ Reads fitted AMISR data from a SRI format HDF5 file

    **Args**:
      * **filepath** (str): The file to be read.

    **Output**: A dictionary with the following keys (where B = beams, N = # time steps,
                and R = # range gates):
                  'az' - B length vector of azimuthal angles in degrees for
                         each ISR beam
                  'el' - B length vector of elevation angles in degrees
                         for each ISR beam
                  'site_latitude' - the ISR site geographic latitude
                  'site_longitude' - the ISR site geographic longitude
                  'site_altitude' - the ISR site altitude above sea level
                                    in meters
                  'times' - Nx2 array of datetime describing the start and
                            end times for an ISR measurement
                  'ave_times' - N array of datetime describing the average
                               time for each ISR measurement
                  'altitude' - BxR array of heights above sea level in
                               meters for each range cell on each beam
                  'range' - BxR array of range along beam direction in
                            meters for each range cell on each beam
                  'density' - NxBxR electron density (not log10)
                  'edensity' - NxBxR error in electron density (not log10)
                  'density_uncor' - NxBxR uncorrected electron density (not log10)
                  'edensity_uncor' - NxBxR error in uncorrected electron density (not log10)
                  'Te' - NxBxR electron temperature in Kelvin
                  'eTe' - NxBxR error in electron density
                  'Ti' - NxBxR ion temperature in Kelvin
                  'eTi' - NxBxR error in ion temperature
                  'vel' - NxBxR line of sight velocity in m/s
                  'evel' - NxBxR error in line of sight velocity
                  'latitude' - BxR geographic latitude of each range
                                      cell
                  'longitude' - BxR geographic longitude of each
                                       range cell
                  'babs' - BxR magnitude of magnetic field in each
                                  range cell
                  'kvec' - Bx3 k-vector of each beam in local North-
                                  East-Up coordinate system.

    **Example**:
      ::
        import visuamisr
        data = visuamisr.read_data('20160302.001_lp_1min-fitcal.h5')

    written by A. S. Reimer, 2013-07
    """
    # So pylint says this function is too long. I'm happy for now.
    # pylint: disable=too-many-statements

    # Provides a consistent API for working with data in fitted files
    # by returning a dictionary of parameters in the fitted files.
    # This makes plotting easier.
    #
    # This function was originally designed to handle the SRI file
    # format only but has a bit of support for Madrigal format hdf5
    # files.

    data = dict()
    with h5py.File(str(filepath), 'r') as h5file:

        if 'BeamCodes' in h5file: #30 sec integrated file
            bckey = '/BeamCodes'
        elif 'Setup/BeamcodeMap' in h5file: #raw samples file
            bckey = '/Setup/BeamcodeMap'
        elif 'Data/Array Layout/1D Parameters/beamid' in h5file: # 2 minute integrated Madrigal file
            bckey = '/Data/Array Layout/1D Parameters/beamid'
        elif 'Data/Table Layout' in h5file: # old 2 minute integrated file
            bckey = None
        else:
            raise ValueError('{} does not conform to expected format.'.format(filepath))

        # get beam information
        if bckey == '/BeamCodes': # SRI file format
            data['beamcodes'] = np.array(h5file['BeamCodes'][:, 0])
            data['az'] = np.array(h5file['BeamCodes'][:, 1])
            data['el'] = np.array(h5file['BeamCodes'][:, 2])
        elif bckey: # Madrigal case
            data['beamcodes'] = np.unique(h5file[bckey][:, 0])
            data['az'] = np.unique(h5file[bckey][:, 1])
            data['el'] = np.unique(h5file[bckey][:, 2])
        else:
            bckey = '/Data/Table Layout'
            data['az'] = np.unique(h5file[bckey]['azm'])
            data['el'] = np.unique(h5file[bckey]['elm'])
            assert data['az'].size == data['el'].size, 'TODO uniquerow'
            data['beamcodes'] = np.arange(data['az'].size)

        # get site information
        try:
            data['site_latitude'] = h5file['Site']['Latitude'][()]
            data['site_longitude'] = h5file['Site']['Longitude'][()]
            data['site_altitude'] = h5file['Site']['Altitude'][()]
            data['site_code'] = h5file['Site']['Code'][()]
        except KeyError:
            pass

        try:
            data['site_name'] = SITECODES[data['site_code']]
        except KeyError:
            data['site_name'] = ''
            print("Site code not in known list, site name not automatically set. "
                  + "You can manually set the name in the data['site_name'] attribute.")

        # get fitted parameters
        if 'FittedParams' in h5file.keys():
            data['density'] = np.array(h5file['FittedParams']['Ne'])
            data['edensity'] = np.array(h5file['FittedParams']['dNe'])
            temp = np.array(h5file['FittedParams']['Fits'])
            data['Te'] = temp[:, :, :, 1, 1]
            data['Ti'] = temp[:, :, :, 0, 1]
            data['vel'] = temp[:, :, :, 1, 3]
            temp = np.array(h5file['FittedParams']['Errors'])
            data['eTe'] = temp[:, :, :, 1, 1]
            data['eTi'] = temp[:, :, :, 0, 1]
            data['evel'] = temp[:, :, :, 1, 3]
            data['range'] = np.array(h5file['FittedParams']['Range'])
            data['altitude'] = np.array(h5file['FittedParams']['Altitude'])
        else:
            data['altitude'] = np.array(h5file['Geomag']['Altitude'])
            print("No fitted data found, so none read.")

        data['density_uncor'] = np.array(h5file['NeFromPower']['Ne_NoTr'])
        data['edensity_uncor'] = np.array(h5file['NeFromPower']['dNeFrac'])
        data['altitude_uncor'] = np.array(h5file['NeFromPower']['Altitude'])

        data['latitude'] = np.array(h5file['Geomag']['Latitude'])
        data['longitude'] = np.array(h5file['Geomag']['Longitude'])

        data['times'] = np.array([[datetime.utcfromtimestamp(x[0]),
                                   datetime.utcfromtimestamp(x[1])]
                                  for x in h5file['Time']['UnixTime']])
        data['ave_times'] = np.array([datetime.utcfromtimestamp((float(x[0]) +
                                                                 float(x[1]))/2)
                                      for x in h5file['Time']['UnixTime']])

        data['babs'] = np.array(h5file['Geomag']['Babs']).T
        data['kvec'] = np.array(h5file['Geomag']['kvec']).T

    return data


####################################################################################
####################################################################################
####################################################################################

class Analyze():
    """ Provides 3 plotting utilities for visualizing AMISR data.

      *********************
      **Class**: Analyze
      *********************
      Read ISR data dictionary from hdf5 file and plot the data.

      **Functions**:
        * :func:`analyze.plot_polar_beam_pattern`
        * :func:`analyze.rti`
        * :func:`analyze.add_title`
        * :func:`analyze.plot_beams3d`
        * :func:`analyze.profile_plot`
    """
    # I'm happy with the number of attributes in this class. Should reconsider
    # if it changes again.
    # pylint: disable=too-many-instance-attributes

    def __init__(self, file_path):
        """ Read in the data to be analyzed/plotted/etc. This method expects
            the file to have been created by the gme.isr.fetchData method.

        **Args**:
          * **file_path** (str): path to a datafile as output by fetchData

        **Example**:
          ::
            import visuamisr
            isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')

        written by A. S. Reimer, 2013-07
        modified by A. S. Reimer 2016-05
        """

        file_path = Path(file_path).expanduser()
        file_name = file_path.name

        #Read the file, but first determine whether it was gzipped or not.
        self.data = read_data(file_path)

        #add the data and some file information to the object
        self.stime = self.data['times'][0, 0]
        self.etime = self.data['times'][-1, 1]
        #grab the instrumnet id from the file name
        self.inst_id = self.data['site_code']
        self.file_path = file_path
        self.file_name = file_name
        #Get site location info
        self.site_lat = self.data['site_latitude']
        self.site_lon = self.data['site_longitude']
        self.site_alt = self.data['site_altitude'] / 1000.0
        (self.num_times, self.num_beams, _) = self.data['density_uncor'].shape


####################################################################################
####################################################################################
####################################################################################


    def plot_polar_beam_pattern(self, min_elevation=None):
        """ Plot the beam positions on a polar plot of azimuth and zenith angle

        **Args**:
          * **[min_elevation]** (int or float): the minimum elevation angle to include in the plot

        **Example**:
          ::
            import visuamisr
            isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
            isr.plot_polar_beam_pattern(min_elevation=10)

        written by A. S. Reimer, 2013-07
        """
        # Pylint complains about too many local variables.
        # pylint: disable=too-many-locals

        assert(not min_elevation or isinstance(min_elevation, (int, float))), \
          "min_elevation must be None, int, or float."
        if not min_elevation:
            min_elevation = np.min(self.data['el'])

        if min_elevation < 10:
            min_elevation = 10
        if min_elevation > 90:
            min_elevation = 90

        #used to set the r axis tick markers
        # note that xmax = np.cos(min_elevation*np.pi/180.0)
        elevations = np.arange(90, 0, -10)
        grid_radii = np.cos(elevations * np.pi / 180.0)
        grid_labels = [str(x) for x in elevations]
        grid_labels[0] = ''   # blank out the 90 degree label
        max_radius = np.cos((min_elevation - 10) * np.pi / 180.0)  # some padding added

        #get elevation angles and azimuth angles of each beam
        beam_rs = np.cos(self.data['el'] * np.pi / 180.0)
        thetas = self.data['az'] * np.pi / 180.0

        #create a polar projection figure and set it up in compass mode
        fig = pyplot.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.75], projection='polar')
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location('N')

        if self.data['site_name'] == SITECODES[92]:
            angle = -180
        else:
            angle = 180
        rgrids = ax.set_rgrids(grid_radii, labels=grid_labels, angle=angle)[1]

        #set the r ticks font size to small
        for rgrid in rgrids:
            rgrid.set_size('small')
        # plot dots for the beams
        ax.plot(thetas, beam_rs, '.r', markersize=10)

        # plot the beam numbers on top of the dots for the beams
        bmnum = [i + 1 for i in range(len(beam_rs))]
        for i, bnum in enumerate(bmnum):
            ax.text(thetas[i], beam_rs[i], str(bnum), weight='bold', zorder=2)
        ax.set_rlim([0, max_radius])

        #plot a title
        text = fig.text(0.5, 0.92, 'Beam Pattern', horizontalalignment='center')
        text.set_size('large')
        fig.show()


####################################################################################
####################################################################################
####################################################################################


    def rti(self, params, time_lim=None, ylim=None, clim=None, cmap=None, bmnum=None,
            use_range=None, show=True):
        """ Create a range time intensity plot

        **Args**:
          * **params** (list): list of strings of parameters to plot: 'density','Te','Ti',
                               'velocity'
          * **[time_lim]** (list): list of datetime.datetime corresponding to start and end times
                                   to plot data from
          * **[ylim]** (list): list of int/float corresponding to range/altitude limits to plot
                               data from
          * **[clim]** (list): list of lists containing the colorbar limits for each parameter
                               plotted
          * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
          * **[bmnum]** (int/float): the beam index of the data to be plotted ie) 5 beam mode has
                                     beams 0-4
          * **[use_range]** (bool): True if Range is to be used for the y-axis instead of altitude.

        **Example**:
          ::
            import visuamisr
            from datetime import datetime
            isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
            isr.rti(['density','Te','Ti','velocity'],
                    time_lim=[datetime(2016,3,2,6,0),datetime(2016,3,2,17)],
                    ylim=[100,500],bmnum=10)


        written by A. S. Reimer, 2013-07
        """
        # These could be cleaned up in the future, but for now ignore.
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements

        #Check inputs
        assert(isinstance(params, list)), "params must be a list of strings"
        for param in params:
            assert(isinstance(param, str) and param in ["density", "Te", "Ti", "velocity"]), \
                "valid parameters are density, Te, Ti, and velocity"
        assert(not time_lim or (isinstance(time_lim, list) and len(time_lim) == 2)), \
            "time_lim must be None or a list of a start and an end time in datetime.datetime"
        if time_lim:
            for lim in time_lim:
                assert(isinstance(lim, datetime)), "time_lim list must be of datetime.datetime"
            assert(time_lim[0] < time_lim[1]), "Wrong order. End time comes before start!"
        assert(not ylim or (isinstance(ylim, list) and len(ylim) == 2)), \
            "ylim must be None or a list of a start and an end range/altitude in int/float"
        if ylim:
            for lim in ylim:
                assert(isinstance(lim, (int, float))), "ylim list entries must be int or float"
            assert(ylim[0] < ylim[1]), "Wrong order. End time comes before start!"
        assert(not clim or (isinstance(clim, list))), \
            "clim must be None or a list of a start and an end value in int/float"
        if clim:
            for lim in clim:
                assert(isinstance(lim[0], (int, float))), "clim list entries must be int or float"
                assert(isinstance(lim[1], (int, float))), "clim list entries must be int or float"
                assert(lim[0] < lim[1]), "Starting values must be smaller than ending values."
        assert(not cmap or isinstance(cmap, (mpl.colors.Colormap, str, list))), \
            "%s, %s, %s " % ("cmap must be None, a matplotlib.colors.Colormap",
                             "or a string describing a matplotlib.colors.Colormap",
                             "or a list of colormaps")
        assert(not bmnum or isinstance(bmnum, (int, float))), "bmnum must be None, int, or float"
        assert(not use_range or isinstance(use_range, bool)), "use_range must be None or bool"

        np.seterr(all='ignore')	#turn off numpy warnings

        #Set some defaults
        if not bmnum:
            bmnum = 0

        if not cmap:
            cmap = 'viridis'

        if isinstance(cmap, str):
            cmaps = [cmap] * len(params)
        else:
            cmaps = cmap

        #grab parameters to be used for RTI
        times = self.data["times"]
        ave_times = self.data["ave_times"]

        if use_range:
            rang = self.data["range"]
            ylabel = 'Range (km)'
        else:
            rang = self.data["altitude"]
            ylabel = 'Altitude (km)'

        if not time_lim:
            time_lim = [times[0, 0], times[-1, -1]]

        if not ylim:
            ylim = [0.0, 800.0]

        #find appropriate time indicies such that we only plot data within time_lim
        tinds = np.where((times[:, 0] >= time_lim[0]) & (times[:, 1] <= time_lim[1]))[0]
        times = times[tinds, :]
        ave_times = ave_times[tinds]

        # #Set up x and y "coordinates" for use with pyplot.fill
        # some values in rang might be nan, especially if altitude
        rng_finite_inds = np.where(np.isfinite(rang[bmnum, :]))[0]

        num_x = times.shape[0]
        num_y = rng_finite_inds.size
        temp_y = rang[bmnum, rng_finite_inds]
        temp_y = np.repeat(temp_y[np.newaxis, :], num_x, axis=0)
        temp_y_diff = np.repeat(np.diff(temp_y[0, :])[np.newaxis, :], num_x, axis=0)
        y_diff = np.zeros(temp_y.shape)
        y_diff[:, 0:-1] = temp_y_diff
        y_diff[:, -1] = temp_y_diff[:, -1]

        # Construct the range array for plotting
        y_plot = np.zeros((num_x+1, num_y+1))
        y_plot[0:-1, 0:-1] = temp_y - y_diff/2
        y_plot[0:-1, -1] = temp_y[:, -1] + y_diff[:, -1]/2
        y_plot[-1, :] = y_plot[-2, :]

        # Construct the time array for plotting
        x_plot = np.ndarray((num_x + 1, num_y + 1), dtype=times.dtype)
        x_plot[:num_x, :] = np.repeat(times[:, 0][:, np.newaxis], num_y + 1, axis=1)
        x_plot[num_x, :] = times[num_x - 1, 1]


        #set up a figure for plotting to
        fig = pyplot.figure(figsize=(11, 8.5))

        #add a title
        az = self.data['az'][bmnum]
        el = self.data['el'][bmnum]
        self.add_title(fig, self.stime, self.data['site_name'], beam=bmnum, az=az, el=el, xmax=.85)

        #iterate through the list of parameters and plot each one
        figtop = .85
        figheight = .8 / len(params)
        for i, param in enumerate(params):
            if param == 'density':
                # Detect if input density is log10 yet or not.
                # If not, make it log10 of density
                if self.data['density'].max() < 10**8:
                    parr = self.data['density']
                else:
                    parr = np.log10(self.data['density'])
                clabel = 'Density\nlog10 /m^3)'
            elif param == 'Te':
                parr = self.data['Te']
                clabel = 'Te (K)'
            elif param == 'Ti':
                parr = self.data['Ti']
                clabel = 'Ti (K)'
            elif param == 'velocity':
                parr = self.data['vel']
                clabel = 'Vlos (m/s)'

            cmap = cmaps[i]

            #calculate the positions of the data axis and colorbar axis
            #for the current parameter and then add them to the figure
            pos = [.1, figtop-figheight*(i+1)+.05, .74, figheight-.04]
            cpos = [.86, figtop-figheight*(i+1)+.05, .03, figheight-.04]
            ax = fig.add_axes(pos)
            cax = fig.add_axes(cpos)

            #set the axis tick markers to face outward
            ax.yaxis.set_tick_params(direction='out')
            ax.xaxis.set_tick_params(direction='out')
            ax.yaxis.set_tick_params(direction='out', which='minor')
            ax.xaxis.set_tick_params(direction='out', which='minor')

            #determine the parameter limits
            if not clim:
                if param == 'density':
                    cbar_lim = [9.0, 12.0]
                elif param == 'Te':
                    cbar_lim = [0.0, 3000.0]
                elif param == 'Ti':
                    cbar_lim = [0.0, 2000.0]
                elif param == 'velocity':
                    cbar_lim = [-500.0, 500.0]
            else:
                cbar_lim = clim[i]

            #only add xtick labels if plotting the last parameter
            if not i == len(params) - 1:
                ax.xaxis.set_ticklabels([])
            else:
                #proper formatting for plotting time as hours and minutes from datetime
                ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
                ax.xaxis_date()
                ax.set_xlabel('UT')

            ax.set_xlim(time_lim)
            ax.set_ylim(ylim)
            ax.set_yticks(np.linspace(ylim[0], ylim[1], num=5))
            ax.set_ylabel(ylabel)

            #plot each data point and color them according to the scalar mapping we created
            num_y = y_plot.shape[1] - 1
            param_to_plot = parr[tinds, bmnum, :num_y]
            ax.pcolormesh(x_plot, y_plot/1000.0, param_to_plot, vmin=cbar_lim[0],
                          vmax=cbar_lim[1], cmap=cmap)

            #add a colorbar and label it properly
            cnorm = colors.Normalize(vmin=cbar_lim[0], vmax=cbar_lim[1])
            cbar = mpl.colorbar.ColorbarBase(cax, norm=cnorm, cmap=cmap)
            cbar.set_label(clabel)
            cbar.set_ticks(np.linspace(cbar_lim[0], cbar_lim[1], num=5))

        #finally show the figure
        if show:
            fig.show()

        #turn warnings back on
        np.seterr(all='warn')


####################################################################################
####################################################################################
####################################################################################

    @staticmethod
    def add_title(fig, date, rad, beam=None, az=None, el=None, datetimes=None, xmin=.1,
                  xmax=.86, y=0.9):
        """draws title for an rti plot

        **Args**:
          * **fig**: the figure object to title
          * **date**: the date being plotted as a datetime object
          * **rad**: the name of the radar
          * **[beam]**: the beam number being plotted
          * **[az]**: the azimuth of the beam being plotted
          * **[el]**: the elevation of the beam being plotted
          * **[datetimes]**: an array of datetimes
          * **[xmin]**: minimum x value to plot in page coords
          * **[xmax]**: maximum x value to plot in page coords
          * **[y]**: y value to put title at in page coords
        * **Returns**:
          *Nothing.

        **Example**:
          ::

            from datetime import datetime
            import visuamisr
            from matplotlib import pyplot
            fig = pyplot.figure()
            visuamisr.add_title(fig,datetime(2011,1,1),'PFISR',beam=7)

        Written by A. S. Reimer 2013/07
        Adapted from rtiTitle in DaViTpy written by AJ 20121002
        """
        # This could be cleaned up in the future, but for now ignore.
        # pylint: disable=too-many-arguments

        fig.text(xmin, y, rad, ha='left', weight=550)

        if datetimes:
            clock_time = ''
            for dtime in datetimes:
                clock_time += '{:s} - '.format(dtime.strftime('%H:%M:%S'))
            clock_time = '{:s} UT'.format(clock_time[0:-2])
            title = '{:d}/{:s}/{:d} {:s}'.format(date.day, calendar.month_name[date.month][:3],
                                                 date.year, clock_time)
            fig.text((xmin + xmax) /2., y, title, weight=550, size='large', ha='center')
        else:
            title = '{:d}/{:s}/{:d}'.format(date.day, calendar.month_name[date.month][:3],
                                            date.year)
            fig.text((xmin + xmax) / 2., y, title, weight=550, size='large', ha='center')

        if not beam is None:
            if not az is None and not el is None:
                fig.text(xmax, y, 'Beam: {:d}\nAz: {:.1f} El: {:.1f}'.format(beam, az, el),
                         weight=550, ha='right')
            else:
                fig.text(xmax, y, 'Beam: {:d}'.format(beam), weight=550, ha='right')


####################################################################################
####################################################################################
####################################################################################


    def profile_plot(self, params, time, param_lim=None, bmnum=None, ylim=None, use_range=None):
        """ Create a profile plot

        **Args**:
          * **params** (list): list of strings of parameters to plot: 'density','Te','Ti',
                              'velocity'
          * **time** (datetime.datetime): the time to plot data for
          * **[param_lim]** (list): list of upper and lower values for each parameter
          * **[bmnum]** (int/float): the beam index of the data to be plotted ie) 5 beam mode has
                                     beams 0-4
          * **[ylim]** (list): list of int/float corresponding to range/altitude limits to plot
                               data from
          * **[use_range]** (bool): True to use range for the y-axis instead of altitude.

        **Example**:
          ::
            import visuamisr
            from datetime import datetime
            isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
            isr.profile_plot(['density','Te','Ti','velocity'],
                             datetime(2016,3,2,14,55),bmnum=10,
                             param_lim=[[10**10,10**12],[0,5000],[0,4000],
                                        [-1000,1000]],use_range=True)



        written by A. S. Reimer, 2013-09
        """
        # These could be cleaned up in the future, but for now ignore.
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements

        #Check inputs
        assert(isinstance(params, list)), "params must be a list of strings"
        for param in params:
            assert(isinstance(param, str) and param in ["density", "Te", "Ti", "velocity"]), \
                "valid parameters are density, Te, Ti, and velocity"

        assert(not ylim or (isinstance(ylim, list) and len(ylim) == 2)), \
            "ylim must be None or a list of a start and an end range/altitude in int/float"
        if ylim:
            for lim in ylim:
                assert(isinstance(lim, (int, float))), "ylim list entries must be int or float"
            assert(ylim[0] < ylim[1]), "Wrong order. End time comes before start!"

        assert(not param_lim or (isinstance(param_lim, list) and len(param_lim) == len(params))), \
            "param_lim must be None or a list of lists of a start and end values in int/float"
        if param_lim:
            for lim in param_lim:
                assert(lim[0] < lim[1]), "Starting values must be smaller than ending values."

        assert(not bmnum or isinstance(bmnum, (int, float))), "bmnum must be None, int, or float"

        assert(not use_range or isinstance(use_range, bool)), "use_range must be None or bool"

        np.seterr(all='ignore') #turn off numpy warnings

        #Set some defaults
        if not bmnum:
            bmnum = 0

        #Get the times that data is available and then determine the index
        #for plotting
        times = self.data['times']
        tinds = np.where(np.logical_and(np.array(time) >= times[:, 0],
                                        np.array(time) <= times[:, 1]))[0].tolist()

        #Now only proceed if the time is found
        if not tinds:
            print("Time not found!")
            return

        tinds = tinds[0]

        if use_range:
            rang = self.data["range"]
            ylabel = 'Range (km)'
        else:
            rang = self.data["altitude"]
            ylabel = 'Altitude (km)'

        if not ylim:
            ylim = [0.0, 800.0]

        #set up a figure for plotting to
        fig = pyplot.figure(figsize=(11, 8.5))

        #add a title
        az = self.data['az'][bmnum]
        el = self.data['el'][bmnum]
        self.add_title(fig, self.stime, self.data['site_name'], beam=bmnum, az=az, el=el,
                       datetimes=[times[tinds, 0], times[tinds, 1]], y=0.92)

        #iterate through the list of parameters and plot each one
        figwidth = .75 / len(params)
        for i, param in enumerate(params):
            if param == 'density':
                parr = self.data['density']
                perr = self.data['edensity']
                plabel = 'Density (/m^3)'
            elif param == 'Te':
                parr = self.data['Te']
                perr = self.data['eTe']
                plabel = 'Te (K)'
            elif param == 'Ti':
                parr = self.data['Ti']
                perr = self.data['eTi']
                plabel = 'Ti (K)'
            elif param == 'velocity':
                parr = self.data['vel']
                perr = self.data['evel']
                plabel = 'Vlos (m/s)'

            #calculate the positions of the data axis and colorbar axis
            #for the current parameter and then add them to the figure
            pos = [.1 + (figwidth + 0.025) * (i), 0.15, figwidth - 0.025, 0.75]

            ax = fig.add_axes(pos)

            if param == 'density':
                ax.set_xscale('log')
            #set the axis tick markers to face outward
            ax.yaxis.set_tick_params(direction='out')
            ax.xaxis.set_tick_params(direction='out')
            ax.yaxis.set_tick_params(direction='out', which='minor')
            ax.xaxis.set_tick_params(direction='out', which='minor')

            #determine the parameter limits
            if not param_lim:
                if param == 'density':
                    plim = [10**10, 10**12]
                elif param == 'Te':
                    plim = [0.0, 3000.0]
                elif param == 'Ti':
                    plim = [0.0, 2000.0]
                elif param == 'velocity':
                    plim = [-500.0, 500.0]
            else:
                plim = param_lim[i]

            #only add xtick labels if plotting the last parameter
            if not i == 0:
                ax.yaxis.set_ticklabels([])
            else:
                ax.set_ylabel(ylabel)

            if i == len(params) - 1:
                ax2 = ax.twinx()
                ax2.set_frame_on(True)
                ax2.patch.set_visible(False)
                ax2.set_ylabel('Ground Range (km)')
                ax2.yaxis.set_tick_params(direction='out')
                ax2.yaxis.set_tick_params(direction='out', which='minor')
                ax2lim = []
                if not use_range:
                    ax2lim.append(ylim[0] / np.tan(np.deg2rad(self.data['el'][bmnum])))
                    ax2lim.append(ylim[1] / np.tan(np.deg2rad(self.data['el'][bmnum])))
                else:
                    ax2lim.append(ylim[0] * np.cos(np.deg2rad(self.data['el'][bmnum])))
                    ax2lim.append(ylim[1] * np.cos(np.deg2rad(self.data['el'][bmnum])))
                ax2.set_ylim(ax2lim)

            ax.set_xlim(plim)
            ax.set_ylim(ylim)
            ax.set_xlabel(plabel)
            numx = np.floor(6 - (len(params) - 1)) if np.floor(6 - (len(params) - 1)) >= 4 else 4
            ax.set_xticks(np.linspace(plim[0], plim[1], num=numx))

            #plot the data
            ax.errorbar(parr[tinds, bmnum, :], rang[bmnum, :] / 1000.0,
                        xerr=perr[tinds, bmnum, :])
            ax.scatter(parr[tinds, bmnum, :], rang[bmnum, :] / 1000.0)

        #finally show the figure
        fig.show()

        #turn warnings back on
        np.seterr(all='warn')



####################################################################################
####################################################################################
####################################################################################


    def plot_beams3d(self, param, time, xlim=None, ylim=None, zmax=None, clim=None, cmap=None,
                     sym_size=5):
        """ Make a plot showing ISR data along each beam in 3D.

        **Args**:
          * **param** (str): The parameter to plot: 'density','Te','Ti','velocity'
          * **time** (datetime.datetime): the time to plot data for
          * **[xlim]** (list): list of int/float corresponding to latitude limits to plot data from
          * **[ylim]** (list): list of int/float corresponding to longitude limits to plot data
                               from
          * **[zmax]** (int/float): maximum altiude to plot data for
          * **[clim]** (list): list of lists containing the colorbar limits for each parameter
                               plotted
          * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
          * **[sym_size]** (int/float): see matplotlib.pyplot.scatter documentation (s parameter)

        **Example**:
          ::
            import visuamisr
            from datetime import datetime
            isr = visuamisr.Analyze('20160302.001_lp_1min-fitcal.h5')
            isr.plot_beams3d('density',datetime(2012,11,24,6,40),sym_size=5, clim=[10,12d)


        written by A. S. Reimer, 2013-08
        """

        # These could be cleaned up in the future, but for now ignore.
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements

        #Check inputs
        assert(isinstance(param, str)), "params must be one of density, Te, Ti, velocity."
        assert(isinstance(time, datetime)), "time must be datetime.datetime"
        assert(not xlim or (isinstance(xlim, list) and len(xlim) == 2)), \
            "xlim must be None or a list of a start and an end Latitude"
        if xlim:
            for lim in xlim:
                assert(isinstance(lim, (int, float))), "xlim list entries must be int or float"
            assert(xlim[0] < xlim[1]), "Starting latitude must be smaller than ending latitude."
        assert(not ylim or (isinstance(ylim, list) and len(ylim) == 2)), \
            "ylim must be None or a list of a start and an end Longitude"
        if ylim:
            for lim in ylim:
                assert(isinstance(lim, (int, float))), "ylim list entries must be int or float"
            assert(ylim[0] < ylim[1]), "Starting longitude must be smaller than ending longitude."

        assert(not zmax or isinstance(zmax, (int, float))), \
            "zmax must be None or the maximum altitude to plot."

        assert(not clim or (isinstance(clim, list) and len(clim) == 2)), \
            "clim must be None or a list of a start and an end value in int/float"
        if clim:
            for lim in clim:
                assert(isinstance(lim, (int, float))), "clim list entries must be int or float"
            assert(clim[0] < clim[1]), "Starting values must be smaller than ending values."

        assert(not cmap or isinstance(cmap, (mpl.colors.Colormap, str, list))), \
            "%s, %s, %s " % ("cmap must be None, a matplotlib.colors.Colormap",
                             "or a string describing a matplotlib.colors.Colormap",
                             "or a list of colormaps")

        np.seterr(all='ignore')	#turn off numpy warnings

        #Use the default colormap if necessary
        if not cmap:
            cmap = 'viridis'

        #Get the times that data is available and then determine the index
        #for plotting
        times = self.data['times']
        tinds = np.where(np.logical_and(np.array(time) >= times[:, 0],
                                        np.array(time) <= times[:, 1]))[0].tolist()

        #Now only proceed if the time is found
        if not tinds:
            print("Time not found!")
            return

        tinds = tinds[0]

        #Get parameter to plot
        lats = self.data['latitude']
        lons = self.data['longitude']
        alts = self.data['altitude'] / 1000.0
        if param == 'density':
            # Detect if input density is log10 yet or not.
            # If not, make it log10 of density
            if self.data['density'].max() < 10**8:
                parr = self.data['density']
            else:
                parr = np.log10(self.data['density'])
            clabel = 'Density log10 /m^3)'
        elif param == 'Te':
            parr = self.data['Te']
            clabel = 'Te (K)'
        elif param == 'Ti':
            parr = self.data['Ti']
            clabel = 'Ti (K)'
        elif param == 'velocity':
            parr = self.data['vel']
            clabel = 'Vlos (m/s)'

        #First create a figure with a 3D projection axis
        fig = pyplot.figure()
        ax = fig.add_axes([0.00, 0.05, 0.74, 0.9], projection='3d')
        ax.patch.set_fill(0)
        cax = fig.add_axes([0.80, 0.2, 0.03, 0.6])

        #set the axis tick markers to face outward
        ax.yaxis.set_tick_params(direction='out')
        ax.xaxis.set_tick_params(direction='out')
        ax.yaxis.set_tick_params(direction='out', which='minor')
        ax.xaxis.set_tick_params(direction='out', which='minor')

        #generate a scalar colormapping to map data to cmap
        if not clim:
            cbar_lim = [9, 12]
        else:
            cbar_lim = clim
        cnorm = colors.Normalize(vmin=cbar_lim[0], vmax=cbar_lim[1])
        scalar_map = cmx.ScalarMappable(norm=cnorm, cmap=cmap)

        #plot the location of the radar
        ax.scatter(self.site_lon, self.site_lat, self.site_alt, s=sym_size, c='black')

        #Add a title
        self.add_title(fig, times[tinds, 0], self.data['site_name'],
                       datetimes=times[tinds, :].tolist())

        #Now plot the data along each beam
        (_, num_beams, num_ranges) = parr.shape
        for bind in range(num_beams):
            for rind in range(num_ranges):
                if np.isfinite(parr[tinds, bind, rind]):
                    ax.scatter(lons[bind, rind], lats[bind, rind], alts[bind, rind], s=sym_size,
                               c=[scalar_map.to_rgba(parr[tinds, bind, rind])],
                               edgecolors=scalar_map.to_rgba(parr[tinds, bind, rind]))

        #set X, Y, and Z limits if necessary
        if not xlim:
            xlim = ax.get_xlim()
        if not ylim:
            ylim = ax.get_ylim()
        if not zmax:
            zmax = 800.0

        #Change number of ticks and their spacing
        ax.set_xticks(np.linspace(xlim[0], xlim[1], num=5))
        for tick in ax.get_xticklabels():
            tick.set_horizontalalignment('right')
        ax.set_yticks(np.linspace(ylim[0], ylim[1], num=5))
        for tick in ax.get_yticklabels():
            tick.set_horizontalalignment('left')
        ax.view_init(elev=20, azim=-60)
        ax.set_zlim([0, zmax])

        #Label the axes
        ax.set_xlabel('\nLongitude')
        ax.set_ylabel('\nLatitude')
        ax.set_zlabel('Altitude')

        #add a colorbar and label it properly
        cbar = mpl.colorbar.ColorbarBase(cax, norm=cnorm, cmap=cmap)
        cbar.set_label(clabel)
        cbar.set_ticks(np.linspace(cbar_lim[0], cbar_lim[1], num=5))

        #show the figure
        fig.show()

        #turn warnings back on
        np.seterr(all='warn')
