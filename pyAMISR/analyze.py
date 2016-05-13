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
              'siteLatitude' - the ISR site geographic latitude
              'siteLongitude' - the ISR site geographic longitude
              'siteAltitude' - the ISR site altitude above sea level 
                                in meters
              'times' - Nx2 array of datetime describing the start and 
                        end times for an ISR measurement
              'aveTimes' - N array of datetime describing the average 
                           time for each ISR measurement
              'altitude' - BxR array of heights above sea level in 
                           meters for each range cell on each beam
              'range' - BxR array of range along beam direction in
                        meters for each range cell on each beam
              'density' - NxBxR electron density (not log10)
              'edensity' - NxBxR error in electron density (not log10)
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
  * :func:`analyze.plotPolarBeamPattern`
  * :func:`analyze.rti`
  * :func:`analyze.addTitle`
  * :func:`analyze.plotBeams3D`
  * :func:`analyze.isrBayesVelocity`
  * :func:`analyze.getBeamGridInds`
  * :func:`analyze.calcHorizSlice`
  * :func:`analyze.getGridCellCorners`
  * :func:`analyze.calcRefractiveIndex`
  * :func:`analyze.overlayData`
  * :func:`analyze.profilePlot`

"""


def read_data(filepath):

    import h5py as h5
    import numpy as np
    from datetime import datetime


    data = dict()
    with h5.File(filepath,'r') as f:

        data['beamcodes'] = np.array(f['BeamCodes'][:,0])
        data['az'] = np.array(f['BeamCodes'][:,1])
        data['el'] = np.array(f['BeamCodes'][:,2])

        data['siteLatitude'] = f['Site']['Latitude'].value
        data['siteLongitude'] = f['Site']['Longitude'].value
        data['siteAltitude'] = f['Site']['Altitude'].value
        data['siteCode'] = f['Site']['Code'].value

#Documentation for the "Fits" entry in the HDF5 file:
#'Fitted parameters, Size: Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, collision frequency, LOS speed), Unit: N/A, Kelvin, s^{-1}, m/s, FLAVOR: numpy'
# So it lists ions first, electrons are the last to be listed
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
        data['latitude'] = np.array(f['Geomag']['Latitude'])
        data['longitude'] = np.array(f['Geomag']['Longitude'])

        data['times'] = np.array([[datetime.utcfromtimestamp(x[0]), datetime.utcfromtimestamp(x[1])] for x in f['Time']['UnixTime']])
        data['aveTimes'] = np.array([datetime.utcfromtimestamp((float(x[0]) + float(x[1]))/2) for x in f['Time']['UnixTime']])

        data['babs'] = np.array(f['Geomag']['Babs']).T
        data['kvec'] = np.array(f['Geomag']['kvec']).T

    return data


####################################################################################
####################################################################################
####################################################################################

class analyze(object):

  def __init__(self,filePath):
    """ Read in the data to be analyzed/plotted/etc. This method expects
        the file to have been created by the gme.isr.fetchData method.
    
    **Args**:    
      * **filePath** (str): path to a datafile as output by fetchData

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
       
    written by A. S. Reimer, 2013-07
    modified by A. S. Reimer 2016-05
    """

    import os
    import gzip, cPickle
    import datetime as dt

    fileName = os.path.split(filePath)[-1].split('.')

    #Read the file, but first determine whether it was gzipped or not.
    self.data = read_data(filePath)

    #add the data and some file information to the object
    self.sTime = self.data['times'][0,0]
    self.eTime = self.data['times'][-1,1]
    #grab the instrumnet id from the file name
    self.instId = self.data['siteCode']
    self.filePath = filePath
    self.fileName = fileName
    #Get site location info
    self.siteLat = self.data['siteLatitude']
    self.siteLon = self.data['siteLongitude']
    self.siteAlt = self.data['siteAltitude']/1000.0 
    (self.numTimes, self.numBeams, self.numRanges)=self.data['density'].shape

####################################################################################
####################################################################################
####################################################################################


  def plotPolarBeamPattern(self, maxZen=None):
    """ Plot the beam positions on a polar plot of azimuth and zenith angle

    **Args**:    
      * **[maxZen]** (int or float): the minimum zenith angle to include in the plot 

    written by A. S. Reimer, 2013-07
    """

    from matplotlib import pyplot
    import numpy as np
    import datetime as dt

    assert(not maxZen or isinstance(maxZen,(int,float))),"minZen must be None, int, or float."
    if not maxZen:
      maxZen = 70.0

    #used to set the r axis tick markers
    radii= range(10,int(maxZen+10),10)

    #get zenith angles and azimuth angles of each beam
    rs=90*np.cos(self.data['el']*np.pi/180.0)
    thetas=self.data['az']*np.pi/180.0

    #create a polar projection figure and set it up in compass mode
    fig=pyplot.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.75], projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_rlim([0,radii[-1]])
    t=ax.set_rgrids(radii,angle=180)[1]
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
    t = fig.text(0.5,0.92,'Beam Pattern', horizontalalignment='center')
    t.set_size('large')
    fig.show()


####################################################################################
####################################################################################
####################################################################################


  def rti(self, params, timeLim=None, yLim=None, cLim=None, cmap=None, bmnum=None, rang=None, show=True):
    """ Create a range time intensity plot
    
    **Args**:    
      * **params** (list): list of strings of parameters to plot: 'density','Te','Ti','velocity'
      * **[timeLim]** (list): list of datetime.datetime corresponding to start and end times to plot data from
      * **[yLim]** (list): list of int/float corresponding to range/altitude limits to plot data from
      * **[cLim]** (list): list of lists containing the colorbar limits for each parameter plotted
      * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
      * **[bmnum]** (int/float): the beam index of the data to be plotted ie) 5 beam mode has beams 0-4
      * **[rang]** (bool): True if Range is to be used for the y-axis instead of Altitude.

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.rti(['density','Te','Ti','velocity'],
                timeLim=[datetime(2012,11,24,6,0),datetime(2012,11,24,7)],
                yLim=[100,500],bmnum=33)
    

    written by A. S. Reimer, 2013-07
    """
    from matplotlib import pyplot
    from matplotlib import dates
    from matplotlib import colors
    import matplotlib as mpl
    import matplotlib.cm as cmx
    import numpy as np
    import datetime as dt

    #Check inputs
    assert(isinstance(params,list)),"params must be a list of strings"
    for p in params: assert(isinstance(p,str) and p in ["density","Te","Ti","velocity"]), \
      "valid parameters are density, Te, Ti, and velocity"

    assert(not timeLim or (isinstance(timeLim,list) and len(timeLim) == 2)), \
      "timeLim must be None or a list of a start and an end time in datetime.datetime"
    if timeLim:
      for l in timeLim: assert(isinstance(l,dt.datetime)),"timeLim list entries must be datetime.datetime"
      assert(timeLim[0] < timeLim[1]),"In this program, we prefer to move forward in time."

    assert(not yLim or (isinstance(yLim,list) and len(yLim) == 2)), \
      "yLim must be None or a list of a start and an end range/altitude in int/float"
    if yLim:
      for l in yLim: assert(isinstance(l,(int,float))),"yLim list entries must be int or float"
      assert(yLim[0] < yLim[1]),"Starting range/altitude must be smaller than ending range/altitude."

    assert(not cLim or (isinstance(cLim,list))), \
      "cLim must be None or a list of a start and an end value in int/float"
    if cLim:
      for l in cLim: 
        assert(isinstance(l[0],(int,float))),"cLim list entries must be int or float"
        assert(isinstance(l[1],(int,float))),"cLim list entries must be int or float"
        assert(l[0] < l[1]),"Starting values must be smaller than ending values."

    assert(not cmap or isinstance(cmap,(matplotlib.colors.Colormap, str))), \
      "cmap must be None, a matplotlib.colors.Colormap, or a string describing a matplotlib.colors.Colormap."

    assert(not bmnum or isinstance(bmnum,(int,float))),"bmnum must be None, int, or float"

    assert(not rang or isinstance(rang, bool)),"rang must be None or bool"

    np.seterr(all='ignore')	#turn off numpy warnings

    #Set some defaults
    if not bmnum:
      bmnum=0

    if not cmap:
      cmap='jet' 
   
    #grab parameters to be used for RTI
    t = self.data["times"]		#array of start and end times for each measurement in datetime.datetime
    cT = self.data["aveTimes"]		#array of time in middle of measurement in datetime.datetime

    if rang:
      r = self.data["range"]		#array of central range of each measurement
      yLabel='Range (km)'
    else:
      r = self.data["altitude"]
      yLabel='Altitude (km)'

    if not timeLim:
      timeLim = [t[0,0],t[-1,-1]]

    if not yLim:
      yLim = [0.0, 800.0]

    #find appropriate time indicies such that we only plot data within timeLim
    tInd = np.where(np.logical_and(t[:,0] >= timeLim[0],t[:,1] <= timeLim[1]))[0]
    t=t[tInd,:]
    cT=cT[tInd]

    #Set up x and y "coordinates" for use with pyplot.fill
    lt = len(t)
    lr = len(r.T)
    x = np.ndarray((lr,lt,4),dtype=dt.datetime)
    y = np.zeros((lr,lt,4))


    #grid the coordinates appropriately
    for i in range(lr):
     for j in range(lt):
       if i==0:
         rStep=r[bmnum,i+1]-r[bmnum,i]
         y[i,j,0]=r[bmnum,0]
         y[i,j,1]=r[bmnum,0]+rStep
         y[i,j,2]=y[i,j,1]
         y[i,j,3]=y[i,j,0]
       else:
         rStep=r[bmnum,i]-r[bmnum,i-1]
         y[i,j,0]=r[bmnum,i-1]
         y[i,j,1]=r[bmnum,i-1]+rStep
         y[i,j,2]=y[i,j,1]
         y[i,j,3]=y[i,j,0]

       x[i,j,0]=t[j,0]
       x[i,j,1]=x[i,j,0]
       x[i,j,2]=t[j,1]
       x[i,j,3]=x[i,j,2]


    #set up a figure for plotting to
    fig = pyplot.figure(figsize=(11,8.5))

    #add a title
    self.addTitle(fig,self.sTime,'RISR-C',bmnum)
    
    #iterate through the list of parameters and plot each one
    figtop = .85
    figheight = .8/len(params)
    for p in range(len(params)):
      if (params[p] == 'density'): 
        #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
        pArr = self.data['density'] if self.data['density'].max() < 10**8 else np.log10(self.data['density'])
        cLabel = 'Density\nlog10 /m^3)'
      elif (params[p] == 'Te'): 
        pArr = self.data['Te']
        cLabel = 'Te (K)'
      elif (params[p] == 'Ti'): 
        pArr = self.data['Ti']
        cLabel = 'Ti (K)'
      elif (params[p] == 'velocity'): 
        pArr = self.data['vel']
        cLabel = 'Vlos (m/s)'

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
      if not cLim:
        if (params[p] == 'density'): cl = [9.0,12.0]
        elif (params[p] == 'Te'): cl = [0.0,3000.0]
        elif (params[p] == 'Ti'): cl = [0.0,2000.0]
        elif (params[p] == 'velocity'): cl = [-500.0,500.0]
      else:
        cl = cLim[p]
      #generate a scalar colormapping to map data to cmap
      cNorm = colors.Normalize(vmin=cl[0], vmax=cl[1])
      scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

      #only add xtick labels if plotting the last parameter
      if not p == len(params)-1:
        ax.xaxis.set_ticklabels([])
      else:
        #proper formatting for plotting time as hours and minutes from datetime
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax.xaxis_date()
        ax.set_xlabel('UT')

      ax.set_xlim(timeLim)
      ax.set_ylim(yLim)
      ax.set_yticks(np.linspace(yLim[0],yLim[1],num=5))
      ax.set_ylabel(yLabel)


      #plot little rectangles for each data point and color them according to the scalar mapping we created
      for i in range(lr):
        for j in range(lt):
          if np.isfinite(pArr[tInd[j],bmnum,i]):
            fills = ax.fill(x[i,j,:],y[i,j,:]/1000.0,color=scalarMap.to_rgba(pArr[tInd[j],bmnum,i]))
      #add a colorbar and label it properly
      cbar = mpl.colorbar.ColorbarBase(cax,norm=cNorm,cmap=cmap)
      cbar.set_label(cLabel)
      cbar.set_ticks(np.linspace(cl[0],cl[1],num=5))
      
    #finally show the figure
    if show:
      fig.show()

    #turn warnings back on
    np.seterr(all='warn')


####################################################################################
####################################################################################
####################################################################################


  def addTitle(self, fig, date, rad, beam=None, time=None, xmin=.1,xmax=.86,y=0.9):
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

        import datetime as dt
        import pyAMISR
        from matplotlib import pyplot
        fig = pyplot.figure()
        pyAMISR.addTitle(fig,dt.datetime(2011,1,1),'PFISR',beam=7)
      
  Written by A. S. Reimer 2013/07
  Adapted from rtiTitle in DaViTpy written by AJ 20121002
  """

    import calendar
    
    fig.text(xmin,y,rad,ha='left',weight=550)

    if time:
      timeStr=''
      for t in time:
        timeStr=timeStr+t.strftime('%H:%M:%S')+' - '
      timeStr=timeStr[0:-2]+'UT'
      title=str(date.day)+'/'+calendar.month_name[date.month][:3]+'/'+str(date.year)+' '+timeStr
      fig.text((xmin+xmax)/2.,y,title, \
        weight=550,size='large',ha='center')
    else:
      fig.text((xmin+xmax)/2.,y,str(date.day)+'/'+calendar.month_name[date.month][:3]+'/'+str(date.year), \
        weight=550,size='large',ha='center')

    if type(beam) != type(None):
      fig.text(xmax,y,'Beam '+str(beam),weight=550,ha='right')
    
      
####################################################################################
####################################################################################
####################################################################################


  def plotBeams3D(self, param, time, xLim=None, yLim=None, zMax=None, cLim=None, cmap=None, symSize=5):
    """ Make a plot showing ISR data along each beam in 3D.
    
    **Args**:    
      * **param** (str): The parameter to plot: 'density','Te','Ti','velocity'
      * **time** (datetime.datetime): the time to plot data for
      * **[xLim]** (list): list of int/float corresponding to latitude limits to plot data from
      * **[yLim]** (list): list of int/float corresponding to longitude limits to plot data from
      * **[zMax]** (int/float): maximum altiude to plot data for
      * **[cLim]** (list): list of lists containing the colorbar limits for each parameter plotted
      * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
      * **[symSize]** (int/float): see matplotlib.pyplot.scatter documentation (s parameter)

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.plotBeams3D('density',datetime(2012,11,24,6,40),symSize=5, cLim=[10,12])
    

    written by A. S. Reimer, 2013-08
    """
    import numpy as np
    from matplotlib import pyplot
    from matplotlib import colors
    import matplotlib as mpl
    import matplotlib.cm as cmx
    from mpl_toolkits.mplot3d import Axes3D
    import datetime as dt

    #Check inputs
    assert(isinstance(param,str)),"params must be one of density, Te, Ti, velocity, refracind, refracindX, or refracindO."
    assert(isinstance(time,dt.datetime)),"time must be datetime.datetime"
    assert(not xLim or (isinstance(xLim,list) and len(xLim) == 2)), \
      "xLim must be None or a list of a start and an end Latitude"
    if xLim:
      for l in xLim: assert(isinstance(l,(int,float))),"xLim list entries must be int or float"
      assert(xLim[0] < xLim[1]),"Starting latitude must be smaller than ending latitude."
    assert(not yLim or (isinstance(yLim,list) and len(yLim) == 2)), \
      "yLim must be None or a list of a start and an end Longitude"
    if yLim:
      for l in yLim: assert(isinstance(l,(int,float))),"yLim list entries must be int or float"
      assert(yLim[0] < yLim[1]),"Starting longitude must be smaller than ending longitude."

    assert(not zMax or isinstance(zMax,(int,float))), \
      "zMax must be None or the maximum altitude to plot."

    assert(not cLim or (isinstance(cLim,list) and len(cLim) == 2)), \
      "cLim must be None or a list of a start and an end value in int/float"
    if cLim:
      for l in cLim: assert(isinstance(l,(int,float))),"cLim list entries must be int or float"
      assert(cLim[0] < cLim[1]),"Starting values must be smaller than ending values."

    assert(not cmap or isinstance(cmap,(matplotlib.colors.Colormap, str))), \
      "cmap must be None, a matplotlib.colors.Colormap, or a string describing a matplotlib.colors.Colormap."
    
    np.seterr(all='ignore')	#turn off numpy warnings

    #Use the default colormap if necessary
    if not cmap:
      cmap='jet'

    #Get the times that data is available and then determine the index
    #for plotting
    times = self.data['times']
    tInd = np.where(np.logical_and(np.array(time) >= times[:,0], np.array(time) <= times[:,1]))[0].tolist()

    #Now only proceed if the time is found
    if (len(tInd) > 0):
      tInd=tInd[0]
 
      #Get parameter to plot
      lats = self.data['latitude']
      lons = self.data['longitude']
      alts = self.data['altitude']/1000.0
      if (param == 'density'): 
        #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
        pArr = self.data['density'] if self.data['density'].max() < 10**8 else np.log10(self.data['density'])
        cLabel = 'Density log10 /m^3)'
      elif (param == 'Te'): 
        pArr = self.data['Te']
        cLabel = 'Te (K)'
      elif (param == 'Ti'): 
        pArr = self.data['Ti']
        cLabel = 'Ti (K)'
      elif (param == 'velocity'): 
        pArr = self.data['vel']
        cLabel = 'Vlos (m/s)'
      elif (param =='refracind'):
        pArr = self.data['refracind']
        cLabel = 'Refractive Index'
      elif (param =='refracindX'):
        pArr = self.data['refracindX']
        cLabel = 'Refractive Index'
      elif (param =='refracindO'):
        pArr = self.data['refracindO']
        cLabel = 'Refractive Index'

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
      if not cLim:
        cl=[9,12]
      else:
        cl=cLim
      cNorm = colors.Normalize(vmin=cl[0], vmax=cl[1])
      scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
  
      #plot the location of the radar
      ax.scatter(self.siteLon, self.siteLat, self.siteAlt, s=symSize, color='black') #scalarMap.to_rgba(10.0))
  
      #Add a title
      self.addTitle(fig, self.sTime, 'RISR-C', time=times[tInd,:].tolist())
 
      #Now plot the data along each beam
      (numT,numB,numR) = pArr.shape
      for b in range(numB):
        for r in range(numR):
          if np.isfinite(pArr[tInd,b,r]):
            ax.scatter(lons[b,r], lats[b,r], alts[b,r], s=symSize, color=scalarMap.to_rgba(pArr[tInd,b,r]))

      #set X, Y, and Z limits if necessary
      if not xLim:
        xLim = ax.get_xlim()
      if not yLim:
        yLim = ax.get_ylim()
      if not zMax:
        zMax = 800.0
  
      #Change number of ticks and their spacing 
      ax.set_xticks(np.linspace(xLim[0],xLim[1],num=5))
      for t in ax.get_xticklabels():
        t.set_horizontalalignment('right')
      ax.set_yticks(np.linspace(yLim[0],yLim[1],num=5))
      for t in ax.get_yticklabels():
        t.set_horizontalalignment('left')
      ax.view_init(elev=20, azim=-60)
      ax.set_zlim([0,zMax])
 
      #Label the axes
      ax.set_xlabel('\nLongitude')
      ax.set_ylabel('\nLatitude')
      ax.set_zlabel('Altitude')
 
      #add a colorbar and label it properly
      cbar = mpl.colorbar.ColorbarBase(cax,norm=cNorm,cmap=cmap)
      cbar.set_label(cLabel)
      cbar.set_ticks(np.linspace(cl[0],cl[1],num=5))
 
      #show the figure
      fig.show()

    else:
      print "Time not found!"

    #turn warnings back on
    np.seterr(all='warn')

  
####################################################################################
####################################################################################
####################################################################################


  def isrBayesVelocity(self,vlosIn,kVecs,beamRanges):
    """ Following Heinselman and Nicolls (2008), calculate a velocity vector and covariance for input line of sight velocities
    
    **Args**:    
      * **vlosIn** (np.array): N column array of line of sight velocities
      * **kVecs** (np.array): Nx3 array of k-vectors of each line of sight velocity
      * **beamRanges** (np.array): N column array of range to each vlosIn

    **Example**:
      ::

        from isrAnalysis import *
        isrAnalysis.isrBayesVelocity(vlosIn,kVecs,beamRanges)
    

    .. note:: For more details on the method, see Heinselman and Nicolls, (2008) RADIO SCIENCE, VOL 43, doi:10.1029/2007RS003805

    written by A. S. Reimer, 2013-07
    """
    import numpy as np

    #First convert ensure that vlosIn is a column vector
    nEl = max(t.shape)
    vlosIn = vlosIn.reshape(nEl,1)

    #Detect and remove any NaNs!
    beamsUsed = np.isfinite(vlosIn);             
    vlosIn = vlosIn[np.where(beamsUsed)]
    kVecs = kVecs[:,np.where(beamsUsed)]
    beamRanges = beamRanges[np.where(beamsUsed)]

    #Build the a priori covariance matricies
    vlosErrorSlope = 1/10.0  #1m/s covariance per 10km altitude (but quadratic relation) (per page 8 Heinselman/Nicolls 2008)
    sigmaE = np.diag((vlosErrorSlope*beamRanges)**2)  #covariance for vlosIn error (could this be calculated from vlos error in ISR data?)
    
    velError=[3000,3000,15]           #covariance for bayesVel as per the paper
    sigmaV=np.diag(velError)          #a priori covariance matrix for bayesVel

    #Calculate the Bayesian velocity (equations 12 and 13)
    A=kVecs

    bayesVel=np.dot(sigmaV, A.T, np.linalg.solve(np.dot(A, sigmaV, A.T) +sigmaE, vlosIn))
    bayesCov=np.linalg.inv(np.dot(A.T, np.linalg.solve(sigmaE, A)) + np.linalg.inv(sigmaV))

    self.beamsUsed = beamsUsed          #an array of the beams the velocities used
    self.bayesVel = bayesVel            #the Bayesian derived velocity vector
    self.bayesCov = bayesCov            #the Bayesian derived covariance matrix for bayesVel


####################################################################################
####################################################################################
####################################################################################


  def getBeamGridInds(self,code):
    """ A function to calculate an array of indicies that can be use to plot ISR data in a grid.
    
    **Args**:    
      * **code** (int/float): the code for each beam pattern

    **Example**:
      ::
        import pyAMISR
        isr=pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.getBeamGridInds(0)
        print isr.beamGridInds

    written by A. S. Reimer, 2013-08
    """

    import numpy as np
    if code==0:
      self.beamGridInds=np.array([[ 7,22,20,23,14], \
                                  [24,25,26,27,13], \
                                  [ 6,28,29,30,31], \
                                  [ 5,32,19,33,12], \
                                  [ 4,34,18,35,11], \
                                  [ 3,36,17,37,10], \
                                  [ 2,38,16,39, 9], \
                                  [ 1,40,15,41, 8]])-1


####################################################################################
####################################################################################
####################################################################################


  def calcHorizSlice(self, param, altitude):
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
        'cornerLat' - the corner latitudes of the beam pattern.
        'cornerLon' - the corner longitude of the beam pattern.

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        interps=isr.calcHorizSlice('density',250.0)

    written by A. S. Reimer, 2013-08
    """

    import numpy as np
    from models import aacgm

    #Get the ISR data to use
    lats = self.data['latitude']
    lons = self.data['longitude']
    alts = self.data['altitude']/1000.0
    if (param == 'density'): 
      pArr = self.data['density']
    elif (param == 'Te'): 
      pArr = self.data['Te']
    elif (param == 'Ti'): 
      pArr = self.data['Ti']
    elif (param == 'velocity'): 
      pArr = self.data['vel']
    elif (param == 'refracind'):
      pArr = self.data['refracind']
    elif (param == 'refracindO'): 
      pArr = self.data['refracindO']
    elif (param == 'refracindX'): 
      pArr = self.data['refracindX']

    #Set up some output arrays
    (numT,numB,numR) = pArr.shape
    pOut = np.zeros((numT,numB))
    latOut = np.zeros((numB))
    lonOut = np.zeros((numB))

    #Loop through each beam
    for b in range(numB):
      #first check to see if requested altitude is equal to exisiting altitude
      rInd_eq = np.where(np.array(altitude) == alts[b,:])[0].tolist()
      if len(rInd_eq) == 0:
        #now get the indicies for the altitude above and below the requested altitude
        rInd_p = np.where(alts[b,:]-np.array(altitude) > 0)[0].tolist()
        rInd_m = np.where(np.array(altitude)-alts[b,:] > 0)[0].tolist()

        if (len(rInd_p) > 0 and len(rInd_m) > 0):
          #if they are found then calculate the weighted average of
          rInd_p=rInd_p[0]
          rInd_m=rInd_m[-1]

          alt_p=alts[b,rInd_p]
          alt_m=alts[b,rInd_m]
          dp=alt_p-altitude
          dm=altitude-alt_m
          dt=dp+dm

          #lats, lons, and data
          pOut[:,b]=pArr[:,b,rInd_p]*dm/dt + pArr[:,b,rInd_m]*dp/dt
          latOut[b]=lats[b,rInd_p]*dm/dt + lats[b,rInd_m]*dp/dt
          lonOut[b]=lons[b,rInd_p]*dm/dt + lons[b,rInd_m]*dp/dt
        else:
          #if no data found, set things to NaN
          pOut[:,b]=np.zeros(numT)*np.nan
          latOut[b]=np.nan
          lonOut[b]=np.nan
      else:
        rInd_eq=rInd_eq[0]
        pOut[:,b]=pArr[:,b,rInd_eq]
        latOut[b]=lats[b,rInd_eq]
        lonOut[b]=lons[b,rInd_eq]

    #Now calculate the corner latitudes and longitude for each beam
    cellCoords = self.getGridCellCorners(latOut,lonOut)
    #Now build a dictionary for output
    outs={}
    outs['data'] = pOut
    outs['lats'] = latOut
    outs['lons'] = lonOut
    outs['cornerLat'] = cellCoords['cornerLat']
    outs['cornerLon'] = cellCoords['cornerLon']

    return outs


####################################################################################
####################################################################################
####################################################################################


  def getGridCellCorners(self,lats,lons):
    """ A function to calculate and return the corners of a grid of input latitudes and longitudes. This should usually only be used by the self.calcHorizSlice method.
    
    **Args**:
      * **lats** : the interpolated latitudes from self.calcHorizSlice method.
      * **lons** : the interpolated longitudes from self.calcHorizSlice method.

    **Output**:
      A dictionary with keys:
        'cornerLat' - the corner latitudes of the beam pattern.
        'cornerLon' - the corner longitude of the beam pattern.

    written by A. S. Reimer, 2013-08
    """
    import numpy as np
    (t,o)=self.beamGridInds.shape
    cornerLat=np.zeros((t,o,4))
    cornerLon=np.zeros((t,o,4))
    t,o=t-1,o-1

    lat=lats[self.beamGridInds] #define for readability
    lon=lons[self.beamGridInds]

    #Now generate the points for the grid
    #INSIDE
    for i in range(1,t):
      for j in range(1,o):
        cornerLat[i,j,0]=(lat[i-1,j-1]+lat[i,j-1]+lat[i,j]+lat[i-1,j])/4
        cornerLat[i,j,1]=(lat[i,j-1]+lat[i+1,j-1]+lat[i+1,j]+lat[i,j])/4
        cornerLat[i,j,2]=(lat[i,j]+lat[i+1,j]+lat[i+1,j+1]+lat[i,j+1])/4
        cornerLat[i,j,3]=(lat[i-1,j]+lat[i,j]+lat[i,j+1]+lat[i-1,j+1])/4
        cornerLon[i,j,0]=(lon[i-1,j-1]+lon[i,j-1]+lon[i,j]+lon[i-1,j])/4
        cornerLon[i,j,1]=(lon[i,j-1]+lon[i+1,j-1]+lon[i+1,j]+lon[i,j])/4
        cornerLon[i,j,2]=(lon[i,j]+lon[i+1,j]+lon[i+1,j+1]+lon[i,j+1])/4
        cornerLon[i,j,3]=(lon[i-1,j]+lon[i,j]+lon[i,j+1]+lon[i-1,j+1])/4
 
    #EDGES
    for i in range(1,t):
      cornerLat[i,0,0]=2*cornerLat[i,1,0]-cornerLat[i,1,3]
      cornerLat[i,0,1]=2*cornerLat[i,1,1]-cornerLat[i,1,2]
      cornerLat[i,0,2]=cornerLat[i,1,1]
      cornerLat[i,0,3]=cornerLat[i,1,0]
  
      cornerLon[i,0,0]=2*cornerLon[i,1,0]-cornerLon[i,1,3]
      cornerLon[i,0,1]=2*cornerLon[i,1,1]-cornerLon[i,1,2]
      cornerLon[i,0,2]=cornerLon[i,1,1]
      cornerLon[i,0,3]=cornerLon[i,1,0]
  					
      cornerLat[i,o,3]=2*cornerLat[i,o-1,3]-cornerLat[i,o-1,0]
      cornerLat[i,o,2]=2*cornerLat[i,o-1,2]-cornerLat[i,o-1,1]
      cornerLat[i,o,1]=cornerLat[i,o-1,2]
      cornerLat[i,o,0]=cornerLat[i,o-1,3]
      cornerLon[i,o,0]=cornerLon[i,o-1,3]
      cornerLon[i,o,1]=cornerLon[i,o-1,2]
      cornerLon[i,o,2]=2*cornerLon[i,o-1,2]-cornerLon[i,o-1,1]
      cornerLon[i,o,3]=2*cornerLon[i,o-1,3]-cornerLon[i,o-1,0]
  

    for i in range(1,o):
      cornerLat[0,i,0]=2*cornerLat[1,i,0]-cornerLat[1,i,1]
      cornerLat[0,i,1]=cornerLat[1,i,0]
      cornerLat[0,i,2]=cornerLat[1,i,3]
      cornerLat[0,i,3]=2*cornerLat[1,i,3]-cornerLat[1,i,2]
      cornerLon[0,i,0]=2*cornerLon[1,i,0]-cornerLon[1,i,1]
      cornerLon[0,i,1]=cornerLon[1,i,0]
      cornerLon[0,i,2]=cornerLon[1,i,3]
      cornerLon[0,i,3]=2*cornerLon[1,i,3]-cornerLon[1,i,2]
      cornerLat[t,i,0]=cornerLat[t-1,i,1]
      cornerLat[t,i,1]=2*cornerLat[t-1,i,1]-cornerLat[t-1,i,0]
      cornerLat[t,i,2]=2*cornerLat[t-1,i,2]-cornerLat[t-1,i,3]
      cornerLat[t,i,3]=cornerLat[t-1,i,2]
      cornerLon[t,i,0]=cornerLon[t-1,i,1]
      cornerLon[t,i,1]=2*cornerLon[t-1,i,1]-cornerLon[t-1,i,0]
      cornerLon[t,i,2]=2*cornerLon[t-1,i,2]-cornerLon[t-1,i,3]
      cornerLon[t,i,3]=cornerLon[t-1,i,2]
    #FIRST CORNER
    cornerLat[0,0,0]=2*cornerLat[1,0,0]-cornerLat[1,0,1]
    cornerLat[0,0,1]=cornerLat[1,0,0]
    cornerLat[0,0,2]=cornerLat[1,1,0]
    cornerLat[0,0,3]=cornerLat[0,1,0]	
    cornerLon[0,0,0]=2*cornerLon[1,0,0]-cornerLon[1,0,1]
    cornerLon[0,0,1]=cornerLon[1,0,0]
    cornerLon[0,0,2]=cornerLon[1,1,0]
    cornerLon[0,0,3]=cornerLon[0,1,0]
    #SECOND CORNER
    cornerLat[t,0,0]=cornerLat[t-1,0,1]
    cornerLat[t,0,1]=2*cornerLat[t-1,0,1]-cornerLat[t-1,0,0]
    cornerLat[t,0,2]=cornerLat[t,1,1]
    cornerLat[t,0,3]=cornerLat[t-1,1,1]
    cornerLon[t,0,0]=cornerLon[t-1,0,1]
    cornerLon[t,0,1]=2*cornerLon[t-1,0,1]-cornerLon[t-1,0,0]
    cornerLon[t,0,2]=cornerLon[t,1,1]
    cornerLon[t,0,3]=cornerLon[t-1,1,1]
    #THIRD CORNER
    cornerLat[t,o,0]=cornerLat[t-1,o-1,2]
    cornerLat[t,o,1]=cornerLat[t,o-1,2]
    cornerLat[t,o,2]=2*cornerLat[t-1,o,2]-cornerLat[t-1,o,3]
    cornerLat[t,o,3]=cornerLat[t-1,o,2]
    cornerLon[t,o,0]=cornerLon[t-1,o-1,2]
    cornerLon[t,o,1]=cornerLon[t,o-1,2]
    cornerLon[t,o,2]=2*cornerLon[t-1,o,2]-cornerLon[t-1,o,3]
    cornerLon[t,o,3]=cornerLon[t-1,o,2]
    #FOURTH CORNER
    cornerLat[0,o,0]=cornerLat[0,o-1,3]
    cornerLat[0,o,1]=cornerLat[1,o-1,3]
    cornerLat[0,o,2]=cornerLat[1,o,3]
    cornerLat[0,o,3]=2*cornerLat[1,o,3]-cornerLat[1,o,2]
    cornerLon[0,o,0]=cornerLon[0,o-1,3]
    cornerLon[0,o,1]=cornerLon[1,o-1,3]
    cornerLon[0,o,2]=cornerLon[1,o,3]
    cornerLon[0,o,3]=2*cornerLon[1,o,3]-cornerLon[1,o,2]

    outLat=np.zeros((self.numBeams,4))
    outLon=np.zeros((self.numBeams,4))
    for i in range(0,self.numBeams):
      (x,y)=np.where(self.beamGridInds==i)
      if len(x):
        outLat[i,:]=cornerLat[x,y,:]
        outLon[i,:]=cornerLon[x,y,:]
      else:
        outLat[i,:]=np.nan
        outLon[i,:]=np.nan
    cellCoords={}
    cellCoords['cornerLat']=outLat
    cellCoords['cornerLon']=outLon

    return cellCoords


####################################################################################
####################################################################################
####################################################################################


  def calcRefractiveIndex(self,freqHF):
    """ A function to calculate the refractive index at HF frequencies from ISR density measurements
    
    **Args**:
      * **freqHF** (int/float): The HF frequency transmitted by an HF radar in MHz.

    **Example**:
      ::
        import pyAMISR
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.calcRefractiveIndex(10.5)

    written by A. S. Reimer, 2013-08
    """

    import numpy as np

    #calculate the angular frequency of the input HF frequency
    angFreqHF=2.0*np.pi*freqHF*10.0**6

    #define some physical constants
    electronCharge=1.602*10.0**(-19)
    electronMass=9.109*10.0**(-31)
    eps0=8.854*10.0**(-12)

    #get the required ISR data for the calculation
    density = self.data['density']
    B = self.data['babs']
    #colFreq = self.data['fits'][:,:,:,5,2]

    #set up the arrays to put results of calculations in
    (numT,numB,numR) = density.shape
    nO=np.zeros((numT,numB,numR),dtype=np.complex)
    nX=np.zeros((numT,numB,numR),dtype=np.complex)

    #iterate through and calculate
    for r in range(numR):
      for b in range(numB):
        for t in range(numT):
          X=(density[t,b,r]*electronCharge**2/(electronMass*eps0))/(angFreqHF)**2
          Y=B[b,r]*electronCharge/(electronMass*angFreqHF)
          Z=0/angFreqHF #include collision frequency in here some day maybe?

          #since HF backscatter occurs when k and B are perpendicular
          theta=np.pi/2.0

          nO[t,b,r]=np.sqrt(    1-X/(  1 - np.complex(0,Z) - (0.5*(Y*np.sin(theta))**2/(1-X-np.complex(0,Z))) + \
            np.sqrt(0.25*(Y*np.sin(theta))**4+(Y*np.cos(theta)*(1-X-np.complex(0,Z)))**2)/(1-X-np.complex(0,Z))  )    )
          nX[t,b,r]=np.sqrt(    1-X/(  1- np.complex(0,Z) - (0.5*(Y*np.sin(theta))**2/(1-X-np.complex(0,Z))) - \
            np.sqrt(0.25*(Y*np.sin(theta))**4+(Y*np.cos(theta)*(1-X-np.complex(0,Z)))**2)/(1-X-np.complex(0,Z))  )    )

    #calculate the average refractive index
    n=(nO+nX)/2.0    

    self.data['refracindO'] = nO
    self.data['refracindX'] = nX
    self.data['refracind'] = n


####################################################################################
####################################################################################
####################################################################################


  def overlayData(self, param, time, altitude, cLim=None, cmap=None, myMap=None, \
                  zorder=3, alpha=1, show=True, colBar=True, colPad=None, sym=None, grid=False, beams=None):
    """ Overlay ISR data at a particular altitude slice onto a basemap.
    
    **Args**:    
      * **param** (str): The parameter to plot: 'density','Te','Ti','velocity, refracind, refracindX, refracindO'
      * **time** (datetime.datetime): the time to plot data for
      * **altitude** (int/float): int/float corresponding to the altitude slice to plot data at
      * **[cLim]** (list): list of lists containing the colorbar limits for each parameter plotted
      * **[cmap]** (matplotlib.colors.Colormap): a colormap to use for each parameter
      * **[myMap]** (utils.mapObj): a colormap to use for each parameter
      * **[zorder]** (int/float): a matplotlib zorder
      * **[alpha]** (int/float): the transparency (0 invisible, 1 fully visible)
      * **[show]** (boolean): whether or not to show the plot
      * **[colBar]** (boolean): plot a colorbar
      * **[colPad]** (str): None or string corresponding to colorbar padding (default '5%').
      * **[sym]** (list): None or list of sym[0]: symbols to plot instead of rectangles and sym[1]: size of symbol
      * **[grid]** (boolean): None or True to specify whether or not to plot a grid
      * **[beams]** (list): list of beams to plot

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.getBeamGridInds(0)
        isr.overlayData('density', datetime(2012,11,24,6,55), 250.0,
                        cLim=[10,12], zorder=4, beams=[36,1,10,4])

    written by A. S. Reimer, 2013-08
    """

    import numpy as np
    import utils
    from matplotlib import pyplot
    from matplotlib import colors
    import matplotlib as mpl
    import matplotlib.cm as cmx
    import numpy as np
    import datetime as dt

    assert(colPad == None or isinstance(colPad,str)),"colPad must be None or a string describing the colorbar padding"
    if not colPad:
      colPad='5%'

    if not beams:
      beams=range(0,self.numBeams)


    #Get the times that data is available and then determine the index
    #for plotting
    times = self.data['times']
    tInd = np.where(np.logical_and(np.array(time) >= times[:,0], np.array(time) <= times[:,1]))[0].tolist()

    #Now only proceed if the time is found
    if (len(tInd) > 0):
      tInd=tInd[0]
      #Get the slice to be plotted
      stuff=self.calcHorizSlice(param, altitude)
      data = stuff['data'][tInd,:]
      cornerLat=stuff['cornerLat'][self.beamGridInds]
      cornerLon=stuff['cornerLon'][self.beamGridInds]


      if (param == 'density'): 
        #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
        pArr = data if data.max() < 10**8 else np.log10(data)
        cLabel = 'Density log10 /m^3)'
      elif (param == 'Te'): 
        pArr = data
        cLabel = 'Te (K)'
      elif (param == 'Ti'): 
        pArr = data
        cLabel = 'Ti (K)'
      elif (param == 'velocity'): 
        pArr = data
        cLabel = 'Vlos (m/s)'
      elif (param in ['refracind','refracindX','refracindO']): 
        pArr = data
        cLabel = 'Refractive Index'
      #determine the parameter limits
      if not cLim:
        if (param == 'density'): cl = [9.0,12.0]
        elif (param == 'Te'): cl = [0.0,3000.0]
        elif (param == 'Ti'): cl = [0.0,2000.0]
        elif (param == 'velocity'): cl = [-500.0,500.0]
        elif (param in ['refracind','refracindX','refracindO']):  cl=[0.5,1.0]
      else:
        cl = cLim

      #determine the color mapping for the data
      if not cmap:
        cmap='jet'
      cNorm = colors.Normalize(vmin=cl[0], vmax=cl[1])
      scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

      #Now we can plot the data on a map
      #if a map object was not passed to this function, create one
      if not myMap:
        fig = pyplot.figure()
        if colBar:
          ax = fig.add_axes([0.05, 0.1, 0.8, 0.8])
        else:
          ax = fig.add_axes([0.1,0.1,0.8,0.8])
        myMap = utils.mapObj(lat_0=self.siteLat,lon_0=self.siteLon,width=1.0e6,height=1.0e6, coords='geo', ax=ax)
      else:
        ax=myMap.ax
        fig=ax.figure

      #plot little rectangles or symbols for each data point and color them according to the scalar mapping we created
      #Symbol stuff
      if sym:
        marker=sym[0]
        if len(sym) > 1: 
          symSize=sym[1]
        else:
          s=20

      #plotting stuff
      bm=self.beamGridInds
      (l,w)=bm.shape
      for i in range(l):
        for j in range(w):
          if np.isfinite(pArr[bm[i,j]]) and bm[i,j] in beams:
            try:
              X,Y=myMap(cornerLon[i,j,:],cornerLat[i,j,:],coords='geo')
            except:
              X,Y=myMap(cornerLon[i,j,:],cornerLat[i,j,:])

            if not sym and not grid:
              fills = ax.fill(X,Y,color=scalarMap.to_rgba(pArr[bm[i,j]]), \
                              zorder=zorder, alpha=alpha, edgecolor='none')
            elif sym:
              ax.scatter(np.average(X),np.average(Y),color=scalarMap.to_rgba(pArr[bm[i,j]]), \
                         marker=marker, s=symSize, zorder=zorder, alpha=alpha, edgecolor=scalarMap.to_rgba(pArr[bm[i,j]]))
            elif grid:
              X=X.tolist()
              Y=Y.tolist()
              X.append(X[0])
              Y.append(Y[0])
              ax.plot(X,Y,'k', zorder=zorder, alpha=alpha)

      #add a colorbar and label it properly
      if colBar:
        if type(colBar) == bool:
          cax, _ = mpl.colorbar.make_axes(ax,location='right')
        else:
          cbar = mpl.colorbar.ColorbarBase(colBar,norm=cNorm,cmap=cmap)
          cbar.set_label(cLabel)
          cbar.set_ticks(np.linspace(cl[0],cl[1],num=5))
        
      #only show the figure if requested (useful for producing many plots)
      if show:
        fig.show()


####################################################################################
####################################################################################
####################################################################################


  def profilePlot(self, params, time, paramLim=None, bmnum=None, yLim=None, rang=None):
    """ Create a profile plot
    
    **Args**:    
      * **params** (list): list of strings of parameters to plot: 'density','Te','Ti','velocity'
      * **time** (datetime.datetime): the time to plot data for
      * **[paramLim]** (list): list of datetime.datetime corresponding to start and end times to plot data from
      * **[bmnum]** (int/float): the beam index of the data to be plotted ie) 5 beam mode has beams 0-4
      * **[yLim]** (list): list of int/float corresponding to range/altitude limits to plot data from
      * **[rang]** (bool): True if Range is to be used for the y-axis instead of Altitude.

    **Example**:
      ::
        import pyAMISR
        from datetime import datetime
        isr = pyAMISR.analyze('20160302.001_lp_1min.h5')
        isr.profilePlot(['density','Te','Ti','velocity'],
                        datetime(2012,11,24,6,5,0),bmnum=40,
                        paramLim=[[10**10,10**12],[0,5000],[0,4000],
                                  [-1000,1000]],rang=True)
    

    written by A. S. Reimer, 2013-09
    """
    from matplotlib import pyplot
    from matplotlib import dates
    from matplotlib import colors
    import matplotlib as mpl
    import matplotlib.cm as cmx
    import numpy as np
    import datetime as dt

    #Check inputs
    assert(isinstance(params,list)),"params must be a list of strings"
    for p in params: assert(isinstance(p,str) and p in ["density","Te","Ti","velocity"]), \
      "valid parameters are density, Te, Ti, and velocity"

    assert(not yLim or (isinstance(yLim,list) and len(yLim) == 2)), \
      "yLim must be None or a list of a start and an end range/altitude in int/float"
    if yLim:
      for l in yLim: assert(isinstance(l,(int,float))),"yLim list entries must be int or float"
      assert(yLim[0] < yLim[1]),"Starting range/altitude must be smaller than ending range/altitude."

    assert(not paramLim or (isinstance(paramLim,list) and len(paramLim) == len(params))), \
      "paramLim must be None or a list of lists of a start and end values in int/float"
    if paramLim:
      for l in paramLim: assert(l[0] < l[1]),"Starting values must be smaller than ending values."

    assert(not bmnum or isinstance(bmnum,(int,float))),"bmnum must be None, int, or float"

    assert(not rang or isinstance(rang, bool)),"rang must be None or bool"

    np.seterr(all='ignore')	#turn off numpy warnings

    #Set some defaults
    if not bmnum:
      bmnum=0

    #Get the times that data is available and then determine the index
    #for plotting
    times = self.data['times']
    tInd = np.where(np.logical_and(np.array(time) >= times[:,0], np.array(time) <= times[:,1]))[0].tolist()

    #Now only proceed if the time is found
    if (len(tInd) > 0):
      tInd=tInd[0]

    #grab parameters to be used for RTI
      t = self.data["times"]		#array of start and end times for each measurement in datetime.datetime
      cT = self.data["aveTimes"]	#array of time in middle of measurement in datetime.datetime

      if rang:
        r = self.data["range"]		#array of central range of each measurement
        yLabel='Range (km)'
      else:
        r = self.data["altitude"]
        yLabel='Altitude (km)'

    if not yLim:
      yLim = [0.0, 800.0]  

    #set up a figure for plotting to
    fig = pyplot.figure(figsize=(11,8.5))

    #add a title
    self.addTitle(fig,self.sTime,'RISR-C',beam=bmnum, time=[times[tInd,0],times[tInd,1]], y=0.92)
    
    #iterate through the list of parameters and plot each one
    figwidth = .75/len(params)
    for p in range(len(params)):
      if (params[p] == 'density'): 
        #detect if input density is log10 yet or not. If not, make it log10 of density (easier to plot)
        pArr = self.data['density']
        pErr = self.data['edensity']
        pLabel = 'Density (/m^3)'
      elif (params[p] == 'Te'): 
        pArr = self.data['Te']
        pErr = self.data['eTe']
        pLabel = 'Te (K)'
      elif (params[p] == 'Ti'): 
        pArr = self.data['Ti']
        pErr = self.data['eTi']
        pLabel = 'Ti (K)'
      elif (params[p] == 'velocity'): 
        pArr = self.data['vel']
        pErr = self.data['evel']
        pLabel = 'Vlos (m/s)'

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
      if not paramLim:
        if (params[p] == 'density'): pLim = [10**10,10**12]
        elif (params[p] == 'Te'): pLim = [0.0,3000.0]
        elif (params[p] == 'Ti'): pLim = [0.0,2000.0]
        elif (params[p] == 'velocity'): pLim = [-500.0,500.0]
      else:
        pLim=paramLim[p]

      #only add xtick labels if plotting the last parameter
      if not p == 0:
        ax.yaxis.set_ticklabels([])
      else:
        ax.set_ylabel(yLabel)

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
          ax2Lim.append(yLim[0]/np.tan(np.deg2rad(self.data['el'][bmnum])))
          ax2Lim.append(yLim[1]/np.tan(np.deg2rad(self.data['el'][bmnum])))
        else:
          ax2Lim.append(yLim[0]*np.cos(np.deg2rad(self.data['el'][bmnum])))
          ax2Lim.append(yLim[1]*np.cos(np.deg2rad(self.data['el'][bmnum])))
        ax2.set_ylim(ax2Lim)



      ax.set_xlim(pLim)
      ax.set_ylim(yLim)
      ax.set_xlabel(pLabel)
      numx = np.floor( 6-(len(params)-1)) if np.floor( 6-(len(params)-1)) >=4 else 4
      ax.set_xticks(np.linspace(pLim[0],pLim[1],num=numx ))

      #plot the data
      fills = ax.errorbar(pArr[tInd,bmnum,:],r[bmnum,:]/1000.0,xerr=pErr[tInd,bmnum,:])
      fills = ax.scatter(pArr[tInd,bmnum,:],r[bmnum,:]/1000.0)
     
    #finally show the figure
    fig.show()

    #turn warnings back on
    np.seterr(all='warn')