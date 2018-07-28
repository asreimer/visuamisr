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


*********************
**Module**: analyze
*********************
Read ISR data dictionary from hdf5 file and plot the data.

**Functions**:
  * :func:`read_data`
  * :func:`analyze.plot_polar_beam_pattern`
  * :func:`analyze.rti`
  * :func:`analyze.add_title`
  * :func:`analyze.plot_beams3D`
  * :func:`analyze.isr_bayes_velocity`
  * :func:`analyze.get_beam_grid_inds`
  * :func:`analyze.calc_horiz_slice`
  * :func:`analyze.get_grid_cell_corners`
  * :func:`analyze.calc_refractive_index`
  * :func:`analyze.profile_plot`

"""
from . import Path
import h5py
import numpy as np
from datetime import datetime
import calendar
from matplotlib import pyplot
from matplotlib import colors
from matplotlib import dates
import matplotlib as mpl
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
#
try:
    import cartopy
except ImportError:
    cartopy = None


# Dictionary of site code number/radar name pairs
site_codes = {91: 'RISR-N', 92: 'RISR-C', 61: 'PFISR'}

def read_data(filepath):

    data = dict()
    with h5py.File(str(filepath),'r') as f:

        if 'BeamCodes' in f: #30 sec integrated file
            bckey = '/BeamCodes'
        elif 'Setup/BeamcodeMap' in f: #raw samples file
            bckey = '/Setup/BeamcodeMap'
        elif 'Data/Array Layout/1D Parameters/beamid' in f: # 2 minute integrated Madrigal file
            bckey = '/Data/Array Layout/1D Parameters/beamid'
        elif 'Data/Table Layout' in f: # old 2 minute integrated file
            bckey = None
        else:
            raise ValueError('{} does not conform to expected format.'.format(filepath))

        # SRI file format        
        if bckey == '/BeamCodes':
            data['beamcodes'] = np.array(f['BeamCodes'][:,0])
            data['az'] = np.array(f['BeamCodes'][:,1])
            data['el'] = np.array(f['BeamCodes'][:,2])

        # unique for Madrigal case
        # Why are we using np.unique here? It breaks things.
        elif bckey: # newer file
            data['beamcodes'] = np.unique(f[bckey][:,0])
            data['az']        = np.unique(f[bckey][:,1])
            data['el']        = np.unique(f[bckey][:,2])
        else: # FIXME old file
            bckey = '/Data/Table Layout'
            data['az']        = np.unique(f[bckey]['azm'])
            data['el']        = np.unique(f[bckey]['elm'])
            assert data['az'].size == data['el'].size,'TODO uniquerow'
            data['beamcodes'] = np.arange(data['az'].size)

        try:
            data['site_latitude'] = f['Site']['Latitude'].value
            data['site_longitude'] = f['Site']['Longitude'].value
            data['site_altitude'] = f['Site']['Altitude'].value
            data['site_code'] = f['Site']['Code'].value
        except KeyError: #FIXME using integrated Madrigal files, it's under /MetaData/Experiment Parameters
            pass

        try:
            data['site_name'] = site_codes[data['site_code']]
        except KeyError:
            data['site_name'] = ''
            print("Site code not in known list, site name not automatically set. "
                  + "You can manually set the name in the data['site_name'] attribute.")

#Documentation for the "Fits" entry in the HDF5 file:
#'Fitted parameters, Size: Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, collision frequency, LOS speed), Unit: N/A, Kelvin, s^{-1}, m/s, FLAVOR: numpy'
# So it lists ions first, electrons are the last to be listed
        if 'FittedParams' in f.keys():
          data['density'] = np.array(f['FittedParams']['Ne'])
          data['edensity'] = np.array(f['FittedParams']['dNe'])
          temp = np.array(f['FittedParams']['Fits'])
          data['Te'] = temp[:,:,:,1,1]
          data['Ti'] = temp[:,:,:,0,1]
          data['vel'] = temp[:,:,:,1,3]
          temp = np.array(f['FittedParams']['Errors'])
          data['eTe'] = temp[:,:,:,1,1]
          data['eTi'] = temp[:,:,:,0,1]
          data['evel'] = temp[:,:,:,1,3]
          data['range'] = np.array(f['FittedParams']['Range'])
          data['altitude'] = np.array(f['FittedParams']['Altitude'])
        else:
          data['altitude'] = np.array(f['Geomag']['Altitude'])
          print("No fitted data found, so none read.")

        data['density_uncor'] = np.array(f['NeFromPower']['Ne_NoTr'])
        data['edensity_uncor'] = np.array(f['NeFromPower']['dNeFrac'])
        data['altitude_uncor'] = np.array(f['NeFromPower']['Altitude'])

        data['latitude'] = np.array(f['Geomag']['Latitude'])
        data['longitude'] = np.array(f['Geomag']['Longitude'])

        data['times'] = np.array([[datetime.utcfromtimestamp(x[0]), datetime.utcfromtimestamp(x[1])] for x in f['Time']['UnixTime']])
        data['ave_times'] = np.array([datetime.utcfromtimestamp((float(x[0]) + float(x[1]))/2) for x in f['Time']['UnixTime']])

        data['babs'] = np.array(f['Geomag']['Babs']).T
        data['kvec'] = np.array(f['Geomag']['kvec']).T

    return data


####################################################################################
####################################################################################
####################################################################################

class analyze(object):

  def __init__(self,file_path):
    """ Read in the data to be analyzed/plotted/etc. This method expects
        the file to have been created by the gme.isr.fetchData method.

    **Args**:
      * **file_path** (str): path to a datafile as output by fetchData

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')

    written by A. S. Reimer, 2013-07
    modified by A. S. Reimer 2016-05
    """

    file_path = Path(file_path).expanduser()

    file_name = file_path.name

    #Read the file, but first determine whether it was gzipped or not.
    self.data = read_data(file_path)

    #add the data and some file information to the object
    self.stime = self.data['times'][0,0]
    self.etime = self.data['times'][-1,1]
    #grab the instrumnet id from the file name
    self.inst_id = self.data['site_code']
    self.file_path = file_path
    self.file_name = file_name
    #Get site location info
    self.site_lat = self.data['site_latitude']
    self.site_lon = self.data['site_longitude']
    self.site_alt = self.data['site_altitude']/1000.0
    (self.num_times, self.num_beams, _)=self.data['density_uncor'].shape

    if not cartopy is None:
      self.default_projection = cartopy.crs.Mercator(central_longitude=self.site_lon,
                                                     min_latitude=self.site_lat-10,
                                                     max_latitude=self.site_lat+10,
                                                     latitude_true_scale=self.site_lat
                                                    )

####################################################################################
####################################################################################
####################################################################################


  def plot_polar_beam_pattern(self, max_zenith=None):
    """ Plot the beam positions on a polar plot of azimuth and zenith angle

    **Args**:
      * **[max_zenith]** (int or float): the minimum zenith angle to include in the plot

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.plot_polar_beam_pattern(max_zenith=85)

    written by A. S. Reimer, 2013-07
    """
    assert(not max_zenith or isinstance(max_zenith,(int,float))),"max_zenith must be None, int, or float."
    if not max_zenith:
      max_zenith = 70.0

    #used to set the r axis tick markers
    radii = range(10,int(max_zenith+10),10)

    #get zenith angles and azimuth angles of each beam
    rs = 90*np.cos(self.data['el']*np.pi/180.0)
    thetas = self.data['az']*np.pi/180.0

    #create a polar projection figure and set it up in compass mode
    fig = pyplot.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.75],projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_rlim([0,radii[-1]])
    t = ax.set_rgrids(radii,angle=180)[1]
    #set the r ticks font size to small
    for i in t:
      i.set_size('small')
    #plot dots for the beams
    ax.plot(thetas,rs,'.r',markersize=10)

    #plot the beam numbers on top of the dots for the beams
    bmnum = [i+1 for i in range(len(rs))]
    for i in range(len(bmnum)):
      ax.text(thetas[i],rs[i],str(bmnum[i]),weight='bold',zorder=2)

    #plot a title
    t = fig.text(0.5,0.92,'Beam Pattern',horizontalalignment='center')
    t.set_size('large')
    fig.show()


####################################################################################
####################################################################################
####################################################################################


  def rti(self, params, time_lim=None, ylim=None, clim=None, cmap=None, bmnum=None, use_range=None, show=True):
    """ Create a range time intensity plot

    **Args**:
      * **params** (list): list of strings of parameters to plot: 'density','Te','Ti','velocity'
      * **[time_lim]** (list): list of datetime.datetime corresponding to start and end times to plot data from
      * **[ylim]** (list): list of int/float corresponding to range/altitude limits to plot data from
      * **[clim]** (list): list of lists containing the colorbar limits for each parameter plotted
      * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
      * **[bmnum]** (int/float): the beam index of the data to be plotted ie) 5 beam mode has beams 0-4
      * **[use_range]** (bool): True if Range is to be used for the y-axis instead of Altitude.

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.rti(['density','Te','Ti','velocity'],
                time_lim=[datetime(2012,11,24,6,0),datetime(2012,11,24,7)],
                ylim=[100,500],bmnum=33)


    written by A. S. Reimer, 2013-07
    """

    #Check inputs
    assert(isinstance(params,list)),"params must be a list of strings"
    for p in params: assert(isinstance(p,str) and p in ["density","Te","Ti","velocity"]), \
      "valid parameters are density, Te, Ti, and velocity"
    assert(not time_lim or (isinstance(time_lim,list) and len(time_lim) == 2)), \
      "time_lim must be None or a list of a start and an end time in datetime.datetime"
    if time_lim:
      for l in time_lim: assert(isinstance(l,datetime)),"time_lim list entries must be datetime.datetime"
      assert(time_lim[0] < time_lim[1]),"In this program, we prefer to move forward in time."
    assert(not ylim or (isinstance(ylim,list) and len(ylim) == 2)), \
      "ylim must be None or a list of a start and an end range/altitude in int/float"
    if ylim:
      for l in ylim: assert(isinstance(l,(int,float))),"ylim list entries must be int or float"
      assert(ylim[0] < ylim[1]),"Starting range/altitude must be smaller than ending range/altitude."
    assert(not clim or (isinstance(clim,list))), \
      "clim must be None or a list of a start and an end value in int/float"
    if clim:
      for l in clim:
        assert(isinstance(l[0],(int,float))),"clim list entries must be int or float"
        assert(isinstance(l[1],(int,float))),"clim list entries must be int or float"
        assert(l[0] < l[1]),"Starting values must be smaller than ending values."
    assert(not cmap or isinstance(cmap,(mpl.colors.Colormap, str, list))), \
      "cmap must be None, a matplotlib.colors.Colormap, or a string describing a matplotlib.colors.Colormap, or a list of colormaps"
    assert(not bmnum or isinstance(bmnum,(int,float))),"bmnum must be None, int, or float"
    assert(not use_range or isinstance(use_range, bool)),"use_range must be None or bool"

    np.seterr(all='ignore')	#turn off numpy warnings

    #Set some defaults
    if not bmnum:
      bmnum = 0

    if not cmap:
      cmap = 'viridis'

    if isinstance(cmap,str):
      cmaps = [cmap] * len(params)
    else:
      cmaps = cmap

    #grab parameters to be used for RTI
    t = self.data["times"]		#array of start and end times for each measurement in datetime.datetime
    ave_times = self.data["ave_times"]		#array of time in middle of measurement in datetime.datetime

    if use_range:
      r = self.data["range"]		#array of central range of each measurement
      ylabel = 'Range (km)'
    else:
      r = self.data["altitude"]
      ylabel = 'Altitude (km)'

    if not time_lim:
      time_lim = [t[0,0],t[-1,-1]]

    if not ylim:
      ylim = [0.0, 800.0]

    #find appropriate time indicies such that we only plot data within time_lim
    tinds = np.where((t[:,0] >= time_lim[0]) & (t[:,1] <= time_lim[1]))[0]
    t = t[tinds,:]
    ave_times = ave_times[tinds]

    #Set up x and y "coordinates" for use with pyplot.fill
    lt = len(t)
    lr = len(r.T)
    x = np.ndarray((lr,lt,4),dtype=datetime)
    y = np.zeros((lr,lt,4))

    #grid the coordinates appropriately
    for i in range(lr):
     for j in range(lt):
       if i==0:
         rstep = r[bmnum,i+1] - r[bmnum,i]
         y[i,j,0] = r[bmnum,0]
         y[i,j,1] = r[bmnum,0] + rstep
         y[i,j,2] = y[i,j,1]
         y[i,j,3] = y[i,j,0]
       else:
         rstep = r[bmnum,i] - r[bmnum,i-1]
         y[i,j,0] = r[bmnum,i-1]
         y[i,j,1] = r[bmnum,i-1] + rstep
         y[i,j,2] = y[i,j,1]
         y[i,j,3] = y[i,j,0]

       x[i,j,0] = t[j,0]
       x[i,j,1] = x[i,j,0]
       x[i,j,2] = t[j,1]
       x[i,j,3] = x[i,j,2]

    #set up a figure for plotting to
    fig = pyplot.figure(figsize=(11,8.5))

    #add a title
    self.add_title(fig,self.stime,self.data['site_name'],bmnum)

    #iterate through the list of parameters and plot each one
    figtop = .85
    figheight = .8 / len(params)
    for p in range(len(params)):
      if (params[p] == 'density'):
        #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
        parr = self.data['density'] if self.data['density'].max() < 10**8 else np.log10(self.data['density'])
        clabel = 'Density\nlog10 /m^3)'
      elif (params[p] == 'Te'):
        parr = self.data['Te']
        clabel = 'Te (K)'
      elif (params[p] == 'Ti'):
        parr = self.data['Ti']
        clabel = 'Ti (K)'
      elif (params[p] == 'velocity'):
        parr = self.data['vel']
        clabel = 'Vlos (m/s)'

      cmap = cmaps[p]

      #calculate the positions of the data axis and colorbar axis
      #for the current parameter and then add them to the figure
      pos = [.1,figtop-figheight*(p+1)+.05,.74,figheight-.04]
      cpos = [.86,figtop-figheight*(p+1)+.05,.03,figheight-.04]
      ax = fig.add_axes(pos)
      cax = fig.add_axes(cpos)

      #set the axis tick markers to face outward
      ax.yaxis.set_tick_params(direction='out')
      ax.xaxis.set_tick_params(direction='out')
      ax.yaxis.set_tick_params(direction='out',which='minor')
      ax.xaxis.set_tick_params(direction='out',which='minor')

      #determine the parameter limits
      if not clim:
        if (params[p] == 'density'): cl = [9.0,12.0]
        elif (params[p] == 'Te'): cl = [0.0,3000.0]
        elif (params[p] == 'Ti'): cl = [0.0,2000.0]
        elif (params[p] == 'velocity'): cl = [-500.0,500.0]
      else:
        cl = clim[p]
      #generate a scalar colormapping to map data to cmap
      cnorm = colors.Normalize(vmin=cl[0],vmax=cl[1])
      scalar_map = cmx.ScalarMappable(norm=cnorm,cmap=cmap)

      #only add xtick labels if plotting the last parameter
      if not (p == len(params) - 1):
        ax.xaxis.set_ticklabels([])
      else:
        #proper formatting for plotting time as hours and minutes from datetime
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax.xaxis_date()
        ax.set_xlabel('UT')

      ax.set_xlim(time_lim)
      ax.set_ylim(ylim)
      ax.set_yticks(np.linspace(ylim[0],ylim[1],num=5))
      ax.set_ylabel(ylabel)

      #plot little rectangles for each data point and color them according to the scalar mapping we created
      for i in range(lr):
        for j in range(lt):
          if np.isfinite(parr[tinds[j],bmnum,i]):
            fills = ax.fill(x[i,j,:],y[i,j,:]/1000.0,color=scalar_map.to_rgba(parr[tinds[j],bmnum,i]))
      #add a colorbar and label it properly
      cbar = mpl.colorbar.ColorbarBase(cax,norm=cnorm,cmap=cmap)
      cbar.set_label(clabel)
      cbar.set_ticks(np.linspace(cl[0],cl[1],num=5))

    #finally show the figure
    if show:
      fig.show()

    #turn warnings back on
    np.seterr(all='warn')


####################################################################################
####################################################################################
####################################################################################


  def add_title(self, fig, date, rad, beam=None, time=None, xmin=.1,xmax=.86,y=0.9):
    """draws title for an rti plot

    **Args**:
      * **fig**: the figure object to title
      * **date**: the date being plotted as a datetime object
      * **rad**: the name of the radar
      * **beam**: the beam number being plotted
      * **[xmin]**: minimum x value to plot in page coords
      * **[xmax]**: maximum x value to plot in page coords
      * **[y]**: y value to put title at in page coords
    * **Returns**:
      *Nothing.

    **Example**:
      ::

        from datetime import datetime
        import pyAMISR
        from matplotlib import pyplot
        fig = pyplot.figure()
        pyAMISR.add_title(fig,datetime(2011,1,1),'PFISR',beam=7)

  Written by A. S. Reimer 2013/07
  Adapted from rtiTitle in DaViTpy written by AJ 20121002
  """
    fig.text(xmin,y,rad,ha='left',weight=550)

    if time:
      timeStr = ''
      for t in time:
        timeStr = timeStr + t.strftime('%H:%M:%S') + ' - '
      timeStr = timeStr[0:-2] + 'UT'
      title = str(date.day) + '/' + calendar.month_name[date.month][:3] + '/' + str(date.year) + ' ' + timeStr
      fig.text((xmin + xmax) /2.,y,title,
               weight=550,size='large',ha='center')
    else:
      fig.text((xmin + xmax) / 2.,y,str(date.day) + '/' + calendar.month_name[date.month][:3] + '/' + str(date.year),
               weight=550,size='large',ha='center')

    if type(beam) != type(None):
      fig.text(xmax,y,'Beam '+ str(beam),weight=550,ha='right')


####################################################################################
####################################################################################
####################################################################################


  def plot_beams3D(self, param, time, xlim=None, ylim=None, zmax=None, clim=None, cmap=None, sym_size=5):
    """ Make a plot showing ISR data along each beam in 3D.

    **Args**:
      * **param** (str): The parameter to plot: 'density','Te','Ti','velocity'
      * **time** (datetime.datetime): the time to plot data for
      * **[xlim]** (list): list of int/float corresponding to latitude limits to plot data from
      * **[ylim]** (list): list of int/float corresponding to longitude limits to plot data from
      * **[zmax]** (int/float): maximum altiude to plot data for
      * **[clim]** (list): list of lists containing the colorbar limits for each parameter plotted
      * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
      * **[sym_size]** (int/float): see matplotlib.pyplot.scatter documentation (s parameter)

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.plot_beams3D('density',datetime(2012,11,24,6,40),sym_size=5, clim=[10,12])


    written by A. S. Reimer, 2013-08
    """

    #Check inputs
    assert(isinstance(param,str)),"params must be one of density, Te, Ti, velocity, refracind, refracindX, or refracindO."
    assert(isinstance(time,datetime)),"time must be datetime.datetime"
    assert(not xlim or (isinstance(xlim,list) and len(xlim) == 2)), \
      "xlim must be None or a list of a start and an end Latitude"
    if xlim:
      for l in xlim: assert(isinstance(l,(int,float))),"xlim list entries must be int or float"
      assert(xlim[0] < xlim[1]),"Starting latitude must be smaller than ending latitude."
    assert(not ylim or (isinstance(ylim,list) and len(ylim) == 2)), \
      "ylim must be None or a list of a start and an end Longitude"
    if ylim:
      for l in ylim: assert(isinstance(l,(int,float))),"ylim list entries must be int or float"
      assert(ylim[0] < ylim[1]),"Starting longitude must be smaller than ending longitude."

    assert(not zmax or isinstance(zmax,(int,float))), \
      "zmax must be None or the maximum altitude to plot."

    assert(not clim or (isinstance(clim,list) and len(clim) == 2)), \
      "clim must be None or a list of a start and an end value in int/float"
    if clim:
      for l in clim: assert(isinstance(l,(int,float))),"clim list entries must be int or float"
      assert(clim[0] < clim[1]),"Starting values must be smaller than ending values."

    assert(not cmap or isinstance(cmap,(mpl.colors.Colormap, str))), \
      "cmap must be None, a matplotlib.colors.Colormap, or a string describing a matplotlib.colors.Colormap."

    np.seterr(all='ignore')	#turn off numpy warnings

    #Use the default colormap if necessary
    if not cmap:
      cmap = 'viridis'

    #Get the times that data is available and then determine the index
    #for plotting
    times = self.data['times']
    tinds = np.where(np.logical_and(np.array(time) >= times[:,0], np.array(time) <= times[:,1]))[0].tolist()

    #Now only proceed if the time is found
    if len(tinds) == 0:
      print("Time not found!")
      return

    tinds = tinds[0]

    #Get parameter to plot
    lats = self.data['latitude']
    lons = self.data['longitude']
    alts = self.data['altitude'] / 1000.0
    if (param == 'density'):
      #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
      parr = self.data['density'] if self.data['density'].max() < 10**8 else np.log10(self.data['density'])
      clabel = 'Density log10 /m^3)'
    elif (param == 'Te'):
      parr = self.data['Te']
      clabel = 'Te (K)'
    elif (param == 'Ti'):
      parr = self.data['Ti']
      clabel = 'Ti (K)'
    elif (param == 'velocity'):
      parr = self.data['vel']
      clabel = 'Vlos (m/s)'
    elif (param =='refracind'):
      parr = self.data['refracind']
      clabel = 'Refractive Index'
    elif (param =='refracindX'):
      parr = self.data['refracindX']
      clabel = 'Refractive Index'
    elif (param =='refracindO'):
      parr = self.data['refracindO']
      clabel = 'Refractive Index'

    #First create a figure with a 3D projection axis
    fig = pyplot.figure()
    ax = fig.add_axes([0.00, 0.05, 0.74, 0.9],projection='3d')
    ax.patch.set_fill(0)
    cax = fig.add_axes([0.80, 0.2,.03, 0.6])

    #set the axis tick markers to face outward
    ax.yaxis.set_tick_params(direction='out')
    ax.xaxis.set_tick_params(direction='out')
    ax.yaxis.set_tick_params(direction='out',which='minor')
    ax.xaxis.set_tick_params(direction='out',which='minor')

    #generate a scalar colormapping to map data to cmap
    if not clim:
      cl = [9,12]
    else:
      cl = clim
    cnorm = colors.Normalize(vmin=cl[0], vmax=cl[1])
    scalar_map = cmx.ScalarMappable(norm=cnorm, cmap=cmap)

    #plot the location of the radar
    ax.scatter(self.site_lon, self.site_lat, self.site_alt, s=sym_size, c='black')

    #Add a title
    self.add_title(fig, times[tinds,0], self.data['site_name'], time=times[tinds,:].tolist())

    #Now plot the data along each beam
    (numT,numB,numR) = parr.shape
    for b in range(numB):
      for r in range(numR):
        if np.isfinite(parr[tinds,b,r]):
          ax.scatter(lons[b,r], lats[b,r], alts[b,r], s=sym_size,
                     c=scalar_map.to_rgba(parr[tinds,b,r]), edgecolors=scalar_map.to_rgba(parr[tinds,b,r]))

    #set X, Y, and Z limits if necessary
    if not xlim:
      xlim = ax.get_xlim()
    if not ylim:
      ylim = ax.get_ylim()
    if not zmax:
      zmax = 800.0

    #Change number of ticks and their spacing
    ax.set_xticks(np.linspace(xlim[0],xlim[1],num=5))
    for t in ax.get_xticklabels():
      t.set_horizontalalignment('right')
    ax.set_yticks(np.linspace(ylim[0],ylim[1],num=5))
    for t in ax.get_yticklabels():
      t.set_horizontalalignment('left')
    ax.view_init(elev=20, azim=-60)
    ax.set_zlim([0,zmax])

    #Label the axes
    ax.set_xlabel('\nLongitude')
    ax.set_ylabel('\nLatitude')
    ax.set_zlabel('Altitude')

    #add a colorbar and label it properly
    cbar = mpl.colorbar.ColorbarBase(cax,norm=cnorm,cmap=cmap)
    cbar.set_label(clabel)
    cbar.set_ticks(np.linspace(cl[0],cl[1],num=5))

    #show the figure
    fig.show()

    #turn warnings back on
    np.seterr(all='warn')


####################################################################################
####################################################################################
####################################################################################


  def isr_bayes_velocity(self,vlos_in,kvecs,beam_ranges):
    """ Following Heinselman and Nicolls (2008), calculate a velocity vector and covariance for input line of sight velocities

    **Args**:
      * **vlos_in** (np.array): N column array of line of sight velocities
      * **kvecs** (np.array): Nx3 array of k-vectors of each line of sight velocity
      * **beam_ranges** (np.array): N column array of range to each vlos_in

    **Example**:
      ::

        from pyAMISR import *
        pyAMISR.isr_bayes_velocity(vlos_in,kvecs,beam_ranges)


    .. note:: For more details on the method, see Heinselman and Nicolls, (2008) RADIO SCIENCE, VOL 43, doi:10.1029/2007RS003805

    written by A. S. Reimer, 2013-07
    """

    #First convert ensure that vlos_in is a column vector
    nEl     = max(t.shape)
    vlos_in = vlos_in.reshape(nEl,1)

    #Detect and remove any NaNs!
    beams_used  = np.isfinite(vlos_in);
    vlos_in     = vlos_in[np.where(beams_used)]
    kvecs       = kvecs[:,np.where(beams_used)]
    beam_ranges = beam_ranges[np.where(beams_used)]

    #Build the a priori covariance matricies
    vlosErrorSlope = 1/10.0  #1m/s covariance per 10km altitude (but quadratic relation) (per page 8 Heinselman/Nicolls 2008)
    sigmaE = np.diag((vlosErrorSlope*beam_ranges)**2)  #covariance for vlos_in error (could this be calculated from vlos error in ISR data?)

    velError = [3000,3000,15]           #covariance for bayes_vel as per the paper
    sigmaV   = np.diag(velError)        #a priori covariance matrix for bayes_vel

    #Calculate the Bayesian velocity (equations 12 and 13)
    A = kvecs

    bayes_vel = np.dot(sigmaV, A.T, np.linalg.solve(np.dot(A, sigmaV, A.T) + sigmaE, vlos_in))
    bayes_cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(sigmaE, A)) + np.linalg.inv(sigmaV))

    self.beams_used = beams_used          #an array of the beams the velocities used
    self.bayes_vel = bayes_vel            #the Bayesian derived velocity vector
    self.bayes_cov = bayes_cov            #the Bayesian derived covariance matrix for bayes_vel


####################################################################################
####################################################################################
####################################################################################


  def get_beam_grid_inds(self,code):
    """ A function to calculate an array of indicies that can be use to plot ISR data in a grid.

    **Args**:
      * **code** (int/float): the code for each beam pattern

    **Example**:
      ::
        import pyAMISR
        isr=pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.get_beam_grid_inds(0)
        print isr.beam_grid_inds

    written by A. S. Reimer, 2013-08
    """

    import numpy as np
    if code == 0:
      #TODO, add more beam patterns
      #Beam Pattern for 40 beam pattern
      grid = np.array([[ 7,22,20,23,14],
                       [24,25,26,27,13],
                       [ 6,28,29,30,31],
                       [ 5,32,19,33,12],
                       [ 4,34,18,35,11],
                       [ 3,36,17,37,10],
                       [ 2,38,16,39, 9],
                       [ 1,40,15,41, 8]]) - 1
    if code == 1:
      #Beam Patter for 42 beam grid within the 51 beam mode
      grid = np.array([[45,44,48,47,46,43,42],
                       [38,37,41,40,39,36,35],
                       [31,30,34,33,32,29,28],
                       [24,23,27,26,25,22,21],
                       [17,16,20,19,18,15,14],
                       [10, 9,13,12,11, 8, 7]]) - 1

    self.beam_grid_inds = grid


####################################################################################
####################################################################################
####################################################################################


  def calc_horiz_slice(self, param, altitude):
    """ A function to calculate and return ISR data, latitude, and longitude at a given altitude. This routine will automatically interpolate the data and coords into a local horizontal plane at the requested altitude.

    **Args**:
      * **param** (str): The parameter to interpolate: 'density', 'Te', 'Ti', 'velocity', 'refracind', 'refracindO', 'refracindX'.
      * **altitude** (int/float): the altitude to calculate the data at in kilometers.
      * **[coords]** (str): either 'geo' or 'mag'

    **Output**:
      A dictionary with keys:
        'data' - the data at requested altitude.
        'lats' - the geographic latitude of the beam.
        'lons' - the geographic longitude of the beam.
        'corner_lat' - the corner latitudes of the beam pattern.
        'corner_lon' - the corner longitude of the beam pattern.

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        interps = isr.calc_horiz_slice('density',250.0)

    written by A. S. Reimer, 2013-08
    """

    import numpy as np

    #Get the ISR data to use
    lats = self.data['latitude']
    lons = self.data['longitude']
    lat_lon_alts = self.data['altitude'] / 1000.0

    if (param == 'density'):
      parr = self.data['density']
    elif (param == 'density_uncor'):
      parr = self.data['density_uncor']
      alts = self.data['altitude_uncor'] / 1000.0
    elif (param == 'Te'):
      parr = self.data['Te']
    elif (param == 'Ti'):
      parr = self.data['Ti']
    elif (param == 'velocity'):
      parr = self.data['vel']
    elif (param == 'refracind'):
      parr = self.data['refracind']
    elif (param == 'refracindO'):
      parr = self.data['refracindO']
    elif (param == 'refracindX'):
      parr = self.data['refracindX']

    #Set up some output arrays
    (numT,numB,numR) = parr.shape
    pout = np.zeros((numT,numB))
    latout = np.zeros((numB))
    lonout = np.zeros((numB))

    #Loop through each beam
    for b in range(numB):
      # Do data first
      #first check to see if requested altitude is equal to exisiting altitude
      rInd_eq = np.where(np.array(altitude) == alts[b,:])[0].tolist()
      if len(rInd_eq) == 0:
        #now get the indicies for the altitude above and below the requested altitude
        rInd_p = np.where(alts[b,:]-np.array(altitude) > 0)[0].tolist()
        rInd_m = np.where(np.array(altitude)-alts[b,:] > 0)[0].tolist()
        if (len(rInd_p) > 0 and len(rInd_m) > 0):
          #if they are found then calculate the weighted average of
          rInd_p = rInd_p[0]
          rInd_m = rInd_m[-1]

          alt_p = alts[b,rInd_p]
          alt_m = alts[b,rInd_m]
          dp = alt_p - altitude
          dm = altitude - alt_m
          dt = dp + dm

          #data
          pout[:,b] = parr[:,b,rInd_p]*dm/dt + parr[:,b,rInd_m]*dp/dt
        else:
          #if no data found, set things to NaN
          pout[:,b] = np.zeros(numT)*np.nan
      else:
        rInd_eq = rInd_eq[0]
        pout[:,b] = parr[:,b,rInd_eq]

      # Now to lats and lons
      #first check to see if requested altitude is equal to exisiting altitude
      rInd_eq = np.where(np.array(altitude) == lat_lon_alts[b,:])[0].tolist()
      if len(rInd_eq) == 0:
        #now get the indicies for the altitude above and below the requested altitude
        rInd_p = np.where(lat_lon_alts[b,:] - np.array(altitude) > 0)[0].tolist()
        rInd_m = np.where(np.array(altitude) - lat_lon_alts[b,:] > 0)[0].tolist()

        if (len(rInd_p) > 0 and len(rInd_m) > 0):
          #if they are found then calculate the weighted average of
          rInd_p = rInd_p[0]
          rInd_m = rInd_m[-1]

          alt_p = lat_lon_alts[b,rInd_p]
          alt_m = lat_lon_alts[b,rInd_m]
          dp = alt_p - altitude
          dm = altitude - alt_m
          dt = dp + dm

          #lats and lons
          latout[b] = lats[b,rInd_p]*dm/dt + lats[b,rInd_m]*dp/dt
          lonout[b] = lons[b,rInd_p]*dm/dt + lons[b,rInd_m]*dp/dt
        else:
          #if no data found, set things to NaN
          latout[b] = np.nan
          lonout[b] = np.nan
      else:
        rInd_eq = rInd_eq[0]
        latout[b] = lats[b,rInd_eq]
        lonout[b] = lons[b,rInd_eq]

    #Now calculate the corner latitudes and longitude for each beam
    cell_coords = self.get_grid_cell_corners(latout,lonout)
    #Now build a dictionary for output
    outs={}
    outs['data'] = pout
    outs['lats'] = latout
    outs['lons'] = lonout
    outs['corner_lat'] = cell_coords['corner_lat']
    outs['corner_lon'] = cell_coords['corner_lon']

    return outs


####################################################################################
####################################################################################
####################################################################################


  def get_grid_cell_corners(self,lats,lons):
    """ A function to calculate and return the corners of a grid of input latitudes and longitudes. This should usually only be used by the self.calc_horiz_slice method.

    **Args**:
      * **lats** : the interpolated latitudes from self.calc_horiz_slice method.
      * **lons** : the interpolated longitudes from self.calc_horiz_slice method.

    **Output**:
      A dictionary with keys:
        'corner_lat' - the corner latitudes of the beam pattern.
        'corner_lon' - the corner longitude of the beam pattern.

    written by A. S. Reimer, 2013-08
    """
    import numpy as np

    (t,o) = self.beam_grid_inds.shape
    corner_lat = np.zeros((t,o,4))
    corner_lon = np.zeros((t,o,4))
    t,o = t-1,o-1

    lat = lats[self.beam_grid_inds] #define for readability
    lon = lons[self.beam_grid_inds]

    #Now generate the points for the grid
    #INSIDE
    for i in range(1,t):
      for j in range(1,o):
        corner_lat[i,j,0] = (lat[i-1,j-1]+lat[i,j-1]+lat[i,j]+lat[i-1,j])/4
        corner_lat[i,j,1] = (lat[i,j-1]+lat[i+1,j-1]+lat[i+1,j]+lat[i,j])/4
        corner_lat[i,j,2] = (lat[i,j]+lat[i+1,j]+lat[i+1,j+1]+lat[i,j+1])/4
        corner_lat[i,j,3] = (lat[i-1,j]+lat[i,j]+lat[i,j+1]+lat[i-1,j+1])/4
        corner_lon[i,j,0] = (lon[i-1,j-1]+lon[i,j-1]+lon[i,j]+lon[i-1,j])/4
        corner_lon[i,j,1] = (lon[i,j-1]+lon[i+1,j-1]+lon[i+1,j]+lon[i,j])/4
        corner_lon[i,j,2] = (lon[i,j]+lon[i+1,j]+lon[i+1,j+1]+lon[i,j+1])/4
        corner_lon[i,j,3] = (lon[i-1,j]+lon[i,j]+lon[i,j+1]+lon[i-1,j+1])/4

    #EDGES
    for i in range(1,t):
      corner_lat[i,0,0] = 2*corner_lat[i,1,0]-corner_lat[i,1,3]
      corner_lat[i,0,1] = 2*corner_lat[i,1,1]-corner_lat[i,1,2]
      corner_lat[i,0,2] = corner_lat[i,1,1]
      corner_lat[i,0,3] = corner_lat[i,1,0]

      corner_lon[i,0,0] = 2*corner_lon[i,1,0]-corner_lon[i,1,3]
      corner_lon[i,0,1] = 2*corner_lon[i,1,1]-corner_lon[i,1,2]
      corner_lon[i,0,2] = corner_lon[i,1,1]
      corner_lon[i,0,3] = corner_lon[i,1,0]

      corner_lat[i,o,3] = 2*corner_lat[i,o-1,3]-corner_lat[i,o-1,0]
      corner_lat[i,o,2] = 2*corner_lat[i,o-1,2]-corner_lat[i,o-1,1]
      corner_lat[i,o,1] = corner_lat[i,o-1,2]
      corner_lat[i,o,0] = corner_lat[i,o-1,3]
      corner_lon[i,o,0] = corner_lon[i,o-1,3]
      corner_lon[i,o,1] = corner_lon[i,o-1,2]
      corner_lon[i,o,2] = 2*corner_lon[i,o-1,2]-corner_lon[i,o-1,1]
      corner_lon[i,o,3] = 2*corner_lon[i,o-1,3]-corner_lon[i,o-1,0]


    for i in range(1,o):
      corner_lat[0,i,0] = 2*corner_lat[1,i,0]-corner_lat[1,i,1]
      corner_lat[0,i,1] = corner_lat[1,i,0]
      corner_lat[0,i,2] = corner_lat[1,i,3]
      corner_lat[0,i,3] = 2*corner_lat[1,i,3]-corner_lat[1,i,2]
      corner_lon[0,i,0] = 2*corner_lon[1,i,0]-corner_lon[1,i,1]
      corner_lon[0,i,1] = corner_lon[1,i,0]
      corner_lon[0,i,2] = corner_lon[1,i,3]
      corner_lon[0,i,3] = 2*corner_lon[1,i,3]-corner_lon[1,i,2]
      corner_lat[t,i,0] = corner_lat[t-1,i,1]
      corner_lat[t,i,1] = 2*corner_lat[t-1,i,1]-corner_lat[t-1,i,0]
      corner_lat[t,i,2] = 2*corner_lat[t-1,i,2]-corner_lat[t-1,i,3]
      corner_lat[t,i,3] = corner_lat[t-1,i,2]
      corner_lon[t,i,0] = corner_lon[t-1,i,1]
      corner_lon[t,i,1] = 2*corner_lon[t-1,i,1]-corner_lon[t-1,i,0]
      corner_lon[t,i,2] = 2*corner_lon[t-1,i,2]-corner_lon[t-1,i,3]
      corner_lon[t,i,3] = corner_lon[t-1,i,2]
    #FIRST CORNER
    corner_lat[0,0,0] = 2*corner_lat[1,0,0]-corner_lat[1,0,1]
    corner_lat[0,0,1] = corner_lat[1,0,0]
    corner_lat[0,0,2] = corner_lat[1,1,0]
    corner_lat[0,0,3] = corner_lat[0,1,0]
    corner_lon[0,0,0] = 2*corner_lon[1,0,0]-corner_lon[1,0,1]
    corner_lon[0,0,1] = corner_lon[1,0,0]
    corner_lon[0,0,2] = corner_lon[1,1,0]
    corner_lon[0,0,3] = corner_lon[0,1,0]
    #SECOND CORNER
    corner_lat[t,0,0] = corner_lat[t-1,0,1]
    corner_lat[t,0,1] = 2*corner_lat[t-1,0,1]-corner_lat[t-1,0,0]
    corner_lat[t,0,2] = corner_lat[t,1,1]
    corner_lat[t,0,3] = corner_lat[t-1,1,1]
    corner_lon[t,0,0] = corner_lon[t-1,0,1]
    corner_lon[t,0,1] = 2*corner_lon[t-1,0,1]-corner_lon[t-1,0,0]
    corner_lon[t,0,2] = corner_lon[t,1,1]
    corner_lon[t,0,3] = corner_lon[t-1,1,1]
    #THIRD CORNER
    corner_lat[t,o,0] = corner_lat[t-1,o-1,2]
    corner_lat[t,o,1] = corner_lat[t,o-1,2]
    corner_lat[t,o,2] = 2*corner_lat[t-1,o,2]-corner_lat[t-1,o,3]
    corner_lat[t,o,3] = corner_lat[t-1,o,2]
    corner_lon[t,o,0] = corner_lon[t-1,o-1,2]
    corner_lon[t,o,1] = corner_lon[t,o-1,2]
    corner_lon[t,o,2] = 2*corner_lon[t-1,o,2]-corner_lon[t-1,o,3]
    corner_lon[t,o,3] = corner_lon[t-1,o,2]
    #FOURTH CORNER
    corner_lat[0,o,0] = corner_lat[0,o-1,3]
    corner_lat[0,o,1] = corner_lat[1,o-1,3]
    corner_lat[0,o,2] = corner_lat[1,o,3]
    corner_lat[0,o,3] = 2*corner_lat[1,o,3]-corner_lat[1,o,2]
    corner_lon[0,o,0] = corner_lon[0,o-1,3]
    corner_lon[0,o,1] = corner_lon[1,o-1,3]
    corner_lon[0,o,2] = corner_lon[1,o,3]
    corner_lon[0,o,3] = 2*corner_lon[1,o,3]-corner_lon[1,o,2]

    outLat = np.zeros((self.num_beams,4))
    outLon = np.zeros((self.num_beams,4))
    for i in range(0,self.num_beams):
      (x,y) = np.where(self.beam_grid_inds == i)
      if len(x):
        outLat[i,:] = corner_lat[x,y,:]
        outLon[i,:] = corner_lon[x,y,:]
      else:
        outLat[i,:] = np.nan
        outLon[i,:] = np.nan
    cell_coords = {}
    cell_coords['corner_lat'] = outLat
    cell_coords['corner_lon'] = outLon

    return cell_coords


####################################################################################
####################################################################################
####################################################################################


  def calc_refractive_index(self,_freq_hf):
    """ A function to calculate the refractive index at HF frequencies from ISR density measurements

    **Args**:
      * **_freq_hf** (int/float): The HF frequency transmitted by an HF radar in MHz.

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.calc_refractive_index(10.5)

    written by A. S. Reimer, 2013-08
    """

    import numpy as np

    #calculate the angular frequency of the input HF frequency
    ang_freq_hf = 2.0*np.pi*freq_hf*10.0**6

    #define some physical constants
    electronCharge = 1.602*10.0**(-19)
    electronMass = 9.109*10.0**(-31)
    eps0 = 8.854*10.0**(-12)

    #get the required ISR data for the calculation
    density = self.data['density']
    B = self.data['babs']
    #colFreq = self.data['fits'][:,:,:,5,2]

    #set up the arrays to put results of calculations in
    (numT,numB,numR) = density.shape
    nO = np.zeros((numT,numB,numR),dtype=np.complex)
    nX = np.zeros((numT,numB,numR),dtype=np.complex)

    #iterate through and calculate
    for r in range(numR):
      for b in range(numB):
        for t in range(numT):
          X = (density[t,b,r]*electronCharge**2/(electronMass*eps0))/(ang_freq_hf)**2
          Y = B[b,r]*electronCharge/(electronMass*ang_freq_hf)
          Z = 0.0/ang_freq_hf #include collision frequency in here some day maybe?

          #since HF backscatter occurs when k and B are perpendicular
          theta = np.pi/2.0

          nO[t,b,r] = np.sqrt(1-X/(1 - np.complex(0,Z) - (0.5*(Y*np.sin(theta))**2/(1-X-np.complex(0,Z))) +
            np.sqrt(0.25*(Y*np.sin(theta))**4+(Y*np.cos(theta)*(1-X-np.complex(0,Z)))**2)/(1-X-np.complex(0,Z))))
          nX[t,b,r] = np.sqrt(1-X/(1- np.complex(0,Z) - (0.5*(Y*np.sin(theta))**2/(1-X-np.complex(0,Z))) -
            np.sqrt(0.25*(Y*np.sin(theta))**4+(Y*np.cos(theta)*(1-X-np.complex(0,Z)))**2)/(1-X-np.complex(0,Z))))

    #calculate the average refractive index
    n = (nO + nX) / 2.0

    self.data['refracindO'] = nO
    self.data['refracindX'] = nX
    self.data['refracind'] = n

####################################################################################
####################################################################################
####################################################################################


  def profile_plot(self, params, time, param_lim=None, bmnum=None, ylim=None, rang=None):
    """ Create a profile plot

    **Args**:
      * **params** (list): list of strings of parameters to plot: 'density','Te','Ti','velocity'
      * **time** (datetime.datetime): the time to plot data for
      * **[param_lim]** (list): list of datetime.datetime corresponding to start and end times to plot data from
      * **[bmnum]** (int/float): the beam index of the data to be plotted ie) 5 beam mode has beams 0-4
      * **[ylim]** (list): list of int/float corresponding to range/altitude limits to plot data from
      * **[rang]** (bool): True if Range is to be used for the y-axis instead of Altitude.

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.profile_plot(['density','Te','Ti','velocity'],
                        datetime(2012,11,24,6,5,0),bmnum=40,
                        param_lim=[[10**10,10**12],[0,5000],[0,4000],
                                   [-1000,1000]],rang=True)


    written by A. S. Reimer, 2013-09
    """

    #Check inputs
    assert(isinstance(params,list)),"params must be a list of strings"
    for p in params: assert(isinstance(p,str) and p in ["density","Te","Ti","velocity"]), \
      "valid parameters are density, Te, Ti, and velocity"

    assert(not ylim or (isinstance(ylim,list) and len(ylim) == 2)), \
      "ylim must be None or a list of a start and an end range/altitude in int/float"
    if ylim:
      for l in ylim: assert(isinstance(l,(int,float))),"ylim list entries must be int or float"
      assert(ylim[0] < ylim[1]),"Starting range/altitude must be smaller than ending range/altitude."

    assert(not param_lim or (isinstance(param_lim,list) and len(param_lim) == len(params))), \
      "param_lim must be None or a list of lists of a start and end values in int/float"
    if param_lim:
      for l in param_lim: assert(l[0] < l[1]),"Starting values must be smaller than ending values."

    assert(not bmnum or isinstance(bmnum,(int,float))),"bmnum must be None, int, or float"

    assert(not rang or isinstance(rang, bool)),"rang must be None or bool"

    np.seterr(all='ignore') #turn off numpy warnings

    #Set some defaults
    if not bmnum:
      bmnum=0

    #Get the times that data is available and then determine the index
    #for plotting
    times = self.data['times']
    tinds = np.where(np.logical_and(np.array(time) >= times[:,0], np.array(time) <= times[:,1]))[0].tolist()

    #Now only proceed if the time is found
    if len(tinds) == 0:
      print("Time not found!")
      return

    tinds=tinds[0]

    #grab parameters to be used for RTI
    t = self.data["times"]    #array of start and end times for each measurement in datetime.datetime
    cT = self.data["ave_times"]  #array of time in middle of measurement in datetime.datetime

    if rang:
      r = self.data["range"]    #array of central range of each measurement
      ylabel = 'Range (km)'
    else:
      r = self.data["altitude"]
      ylabel = 'Altitude (km)'

    if not ylim:
      ylim = [0.0, 800.0]

    #set up a figure for plotting to
    fig = pyplot.figure(figsize=(11,8.5))

    #add a title
    self.add_title(fig,self.stime, self.data['site_name'],beam=bmnum, time=[times[tinds,0],times[tinds,1]], y=0.92)

    #iterate through the list of parameters and plot each one
    figwidth = .75/len(params)
    for p in range(len(params)):
      if (params[p] == 'density'):
        #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
        parr = self.data['density']
        pErr = self.data['edensity']
        plabel = 'Density (/m^3)'
      elif (params[p] == 'Te'):
        parr = self.data['Te']
        pErr = self.data['eTe']
        plabel = 'Te (K)'
      elif (params[p] == 'Ti'):
        parr = self.data['Ti']
        pErr = self.data['eTi']
        plabel = 'Ti (K)'
      elif (params[p] == 'velocity'):
        parr = self.data['vel']
        pErr = self.data['evel']
        plabel = 'Vlos (m/s)'

      #calculate the positions of the data axis and colorbar axis
      #for the current parameter and then add them to the figure
      pos = [.1 + (figwidth+0.025)*(p), 0.15, figwidth-0.025, 0.75]

      ax = fig.add_axes(pos)

      if (params[p] == 'density'): ax.set_xscale('log')
      #set the axis tick markers to face outward
      ax.yaxis.set_tick_params(direction='out')
      ax.xaxis.set_tick_params(direction='out')
      ax.yaxis.set_tick_params(direction='out',which='minor')
      ax.xaxis.set_tick_params(direction='out',which='minor')

      #determine the parameter limits
      if not param_lim:
        if (params[p] == 'density'): plim = [10**10,10**12]
        elif (params[p] == 'Te'): plim = [0.0,3000.0]
        elif (params[p] == 'Ti'): plim = [0.0,2000.0]
        elif (params[p] == 'velocity'): plim = [-500.0,500.0]
      else:
        plim = param_lim[p]

      #only add xtick labels if plotting the last parameter
      if not p == 0:
        ax.yaxis.set_ticklabels([])
      else:
        ax.set_ylabel(ylabel)

      if p == len(params)-1:
        ax2=ax.twinx()
        #ax2.spines['right'].set_position(('axes', 1.2))
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        ax2.set_ylabel('Ground Range (km)')
        ax2.yaxis.set_tick_params(direction='out')
        ax2.yaxis.set_tick_params(direction='out',which='minor')
        ax2Lim=[]
        if not rang:
          ax2Lim.append(ylim[0]/np.tan(np.deg2rad(self.data['el'][bmnum])))
          ax2Lim.append(ylim[1]/np.tan(np.deg2rad(self.data['el'][bmnum])))
        else:
          ax2Lim.append(ylim[0]*np.cos(np.deg2rad(self.data['el'][bmnum])))
          ax2Lim.append(ylim[1]*np.cos(np.deg2rad(self.data['el'][bmnum])))
        ax2.set_ylim(ax2Lim)



      ax.set_xlim(plim)
      ax.set_ylim(ylim)
      ax.set_xlabel(plabel)
      numx = np.floor( 6-(len(params)-1)) if np.floor( 6-(len(params)-1)) >=4 else 4
      ax.set_xticks(np.linspace(plim[0],plim[1],num=numx ))

      #plot the data
      fills = ax.errorbar(parr[tinds,bmnum,:],r[bmnum,:]/1000.0,xerr=pErr[tinds,bmnum,:])
      fills = ax.scatter(parr[tinds,bmnum,:],r[bmnum,:]/1000.0)

    #finally show the figure
    fig.show()

    #turn warnings back on
    np.seterr(all='warn')



####################################################################################
####################################################################################
####################################################################################

# TEMPORARILY REMOVED UNTIL CARTOPY IS INTEGRATED FULLY. REQUIRES A MAJOR REWRITE
  # def overlayData(self, param, time, altitude, clim=None, cmap=None, myax=None,
  #                 zorder=3, alpha=1, show=True, colBar=True, colPad=None, sym=None,
  #                 grid=False, beams=None):
  #   """ Overlay ISR data at a particular altitude slice onto a basemap.

  #   **Args**:
  #     * **param** (str): The parameter to plot: 'density','Te','Ti','velocity, refracind, refracindX, refracindO'
  #     * **time** (datetime.datetime): the time to plot data for
  #     * **altitude** (int/float): int/float corresponding to the altitude slice to plot data at
  #     * **[clim]** (list): list of lists containing the colorbar limits for each parameter plotted
  #     * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
  #     * **[myax]** (matplotlib.figure): a matplotlib figure object
  #     * **[zorder]** (int/float): a matplotlib zorder
  #     * **[alpha]** (int/float): the transparency (0 invisible, 1 fully visible)
  #     * **[show]** (boolean): whether or not to show the plot
  #     * **[colBar]** (boolean): plot a colorbar
  #     * **[colPad]** (str): None or string corresponding to colorbar padding (default '5%').
  #     * **[sym]** (list): None or list of sym[0]: symbols to plot instead of rectangles and sym[1]: size of symbol
  #     * **[grid]** (boolean): None or True to specify whether or not to plot a grid
  #     * **[beams]** (list): list of beams to plot

  #   **Example**:
  #     ::
  #       import pyAMISR
  #       from datetime import datetime
  #       isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
  #       isr.get_beam_grid_inds(0)
  #       isr.overlayData('density', datetime(2012,11,24,6,55), 250.0,
  #                       clim=[10,12], zorder=4, beams=[36,1,10,4])

  #   written by A. S. Reimer, 2013-08
  #   """

  #   assert(colPad == None or isinstance(colPad,str)),"colPad must be None or a string describing the colorbar padding"
  #   if not colPad:
  #     colPad='5%'

  #   if not beams:
  #     beams=range(0,self.num_beams)


  #   #Get the times that data is available and then determine the index
  #   #for plotting
  #   times = self.data['times']
  #   tinds = np.where(np.logical_and(np.array(time) >= times[:,0], np.array(time) <= times[:,1]))[0].tolist()

  #   #Now only proceed if the time is found
  #   if (len(tinds) > 0):
  #     tinds=tinds[0]
  #     #Get the slice to be plotted
  #     stuff=self.calc_horiz_slice(param, altitude)
  #     data = stuff['data'][tinds,:]
  #     corner_lat=stuff['corner_lat'][self.beam_grid_inds]
  #     corner_lon=stuff['corner_lon'][self.beam_grid_inds]


  #     if (param == 'density'):
  #       #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
  #       parr = data if data.max() < 10**8 else np.log10(data)
  #       clabel = 'Density log10 /m^3)'
  #     if (param == 'density_uncor'):
  #       #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
  #       parr = data if data.max() < 10**8 else np.log10(data)
  #       clabel = 'Uncorrected Density log10 /m^3)'
  #     elif (param == 'Te'):
  #       parr = data
  #       clabel = 'Te (K)'
  #     elif (param == 'Ti'):
  #       parr = data
  #       clabel = 'Ti (K)'
  #     elif (param == 'velocity'):
  #       parr = data
  #       clabel = 'Vlos (m/s)'
  #     elif (param in ['refracind','refracindX','refracindO']):
  #       parr = data
  #       clabel = 'Refractive Index'
  #     #determine the parameter limits
  #     if not clim:
  #       if (param == 'density'): cl = [9.0,12.0]
  #       elif (param == 'Te'): cl = [0.0,3000.0]
  #       elif (param == 'Ti'): cl = [0.0,2000.0]
  #       elif (param == 'velocity'): cl = [-500.0,500.0]
  #       elif (param in ['refracind','refracindX','refracindO']):  cl=[0.5,1.0]
  #     else:
  #       cl = clim

  #     #determine the color mapping for the data
  #     if not cmap:
  #       cmap='jet'
  #     cnorm = colors.Normalize(vmin=cl[0], vmax=cl[1])
  #     scalar_map = cmx.ScalarMappable(norm=cnorm, cmap=cmap)

  #     #Now we can plot the data on a map
  #     #if a map object was not passed to this function, create one
  #     if not myax:
  #       fig = pyplot.figure()
  #       ax = fig.add_subplot(111,projection=default_projection)
  #     else:
  #       ax = myax

  #     #plot little rectangles or symbols for each data point and color them according to the scalar mapping we created
  #     #Symbol stuff
  #     if sym:
  #       marker=sym[0]
  #       if len(sym) > 1:
  #         sym_size=sym[1]
  #       else:
  #         s=20

  #     #plotting stuff
  #     bm=self.beam_grid_inds
  #     (l,w)=bm.shape
  #     for i in range(l):
  #       for j in range(w):
  #         if np.isfinite(parr[bm[i,j]]) and bm[i,j] in beams:
  #           try:
  #             X,Y=myMap(corner_lon[i,j,:],corner_lat[i,j,:],coords='geo')
  #           except:
  #             X,Y=myMap(corner_lon[i,j,:],corner_lat[i,j,:])

  #           if not sym and not grid:
  #             fills = ax.fill(X,Y,color=scalar_map.to_rgba(parr[bm[i,j]]), \
  #                             zorder=zorder, alpha=alpha, edgecolor='none')
  #           elif sym:
  #             ax.scatter(np.average(X),np.average(Y),color=scalar_map.to_rgba(parr[bm[i,j]]), \
  #                        marker=marker, s=sym_size, zorder=zorder, alpha=alpha, edgecolor=scalar_map.to_rgba(parr[bm[i,j]]))
  #           elif grid:
  #             X=X.tolist()
  #             Y=Y.tolist()
  #             X.append(X[0])
  #             Y.append(Y[0])
  #             ax.plot(X,Y,'k', zorder=zorder, alpha=alpha)

  #     #add a colorbar and label it properly
  #     if colBar:
  #       if type(colBar) == bool:
  #         cax, _ = mpl.colorbar.make_axes(ax,location='right')
  #       else:
  #         cax = colBar
  #       cbar = mpl.colorbar.ColorbarBase(cax,norm=cnorm,cmap=cmap)
  #       cbar.set_label(clabel)
  #       cbar.set_ticks(np.linspace(cl[0],cl[1],num=5))

  #     #only show the figure if requested (useful for producing many plots)
  #     if show:
  #       fig.show()


####################################################################################
####################################################################################
####################################################################################

# TEMPORARILY REMOVED UNTIL CARTOPY IS INTEGRATED FULLY. REQUIRES A MAJOR REWRITE

  # def overlayBeamGrid(self, altitude, myMap=None, fill=False, fillColor='blue', sym=None, symColor='black', zorder=3, alpha=1):
  #   """ Overlay horizontal beam grid at a particular altitude slice onto a basemap.

  #   **Args**:
  #     * **altitude** (int/float): int/float corresponding to the altitude slice to plot data at
  #     * **[myMap]** (utils.mapObj): a colormap to use for each parameter
  #     * **[fill]** (bool): Specify whether or not to fill the grid with a colour
  #     * **[fillColor]** (str): A string describing the colour that should be used to fill with
  #     * **[sym]** (list): None or list of sym[0]: symbols to plot instead of rectangles and sym[1]: size of symbol
  #     * **[symColor]** (str): A string describing the colour of the symbol to plot
  #     * **[zorder]** (int/float): a matplotlib zorder
  #     * **[alpha]** (int/float): the transparency (0 invisible, 1 fully visible)

  #   **Example**:
  #     ::
  #       import pyAMISR
  #       from datetime import datetime
  #       isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
  #       isr.get_beam_grid_inds(0)
  #       isr.overlayBeamGrid(250.0)

  #   written by A. S. Reimer, 2016-07
  #   """


  #   beams=range(0,self.num_beams)

  #   #Get the slice to be plotted
  #   temp = self.calc_horiz_slice('density', altitude)
  #   corner_lat=temp['corner_lat'][self.beam_grid_inds]
  #   corner_lon=temp['corner_lon'][self.beam_grid_inds]

  #   #Now we can plot the data on a map
  #   #if a map object was not passed to this function, create one
  #   show = False
  #   if not myMap:
  #     fig = pyplot.figure()
  #     ax = fig.add_axes([0.1,0.1,0.8,0.8])
  #     myMap = utils.mapObj(lat_0=self.site_lat,lon_0=self.site_lon,
  #                          width=1.0e6,height=1.0e6,coords='geo',ax=ax)
  #     show = True
  #   else:
  #     ax=myMap.ax
  #     fig=ax.figure

  #   #plot little rectangles or symbols for each data point and color them according to the scalar mapping we created
  #   #Symbol stuff
  #   if sym:
  #     marker=sym[0]
  #     if len(sym) > 1:
  #       sym_size=sym[1]
  #     else:
  #       sym_size=20

  #   #plot the grid
  #   bm=self.beam_grid_inds
  #   (l,w)=bm.shape
  #   for i in range(l):
  #     for j in range(w):
  #       try:
  #         X,Y=myMap(corner_lon[i,j,:],corner_lat[i,j,:],coords='geo')
  #       except:
  #         X,Y=myMap(corner_lon[i,j,:],corner_lat[i,j,:])
  #       if sym:
  #         myMap.scatter(np.average(X),np.average(Y),color=symColor,
  #                     marker=marker, s=sym_size, zorder=zorder,
  #                     alpha=alpha, edgecolor=symColor,latlon=False)
  #       else:
  #         if fill:
  #           fills = ax.fill(X,Y,color=fillColor, \
  #                           zorder=zorder, alpha=alpha, edgecolor='none')
  #         X=X.tolist()
  #         Y=Y.tolist()
  #         X.append(X[0])
  #         Y.append(Y[0])
  #         myMap.plot(X,Y,'k',latlon=False, zorder=zorder, alpha=alpha)
  #   if show:
  #     fig.show()


