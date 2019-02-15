#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 12:01:37 2019

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np

from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import pickle
import glob
import sys,getopt
from matplotlib.colors import LinearSegmentedColormap


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

#file = '20170826-013200.netcdf'

file = '20170823-223509.netcdf'

nc = Dataset(file, 'r')

rotation = nc.variables['RotationTrack30min'][:]
lat = np.linspace(20.0005,55.0005,7000)
lon = np.linspace(-130.0005,-60.0005,14000)
lat2,lon2 = np.meshgrid(lat,lon)





plt.figure(figsize=(18,12))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-99,-96]); ylim = np.array([25,28])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
   
cs = m.contour(lon2,lat2,rotation.T, cmap = 'YlOrRd') 
m.drawcounties()

plt.show()

