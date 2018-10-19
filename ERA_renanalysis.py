#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:43:58 2018

@author: noahbrauer
"""

from ecmwfapi import ECMWFDataServer
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator



import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

file = 'era.nc'
file2 = 'shear.nc'


era = {}

nc = Dataset(file, 'r', unpack = True)
nc2 = Dataset(file2, 'r', unpack = True)


lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:] - 360
lat2, lon2 = np.meshgrid(lat,lon)[:]

time = nc.variables['time'][:]
timeUnits = nc.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
era['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
era['day'] = np.asarray([d.day for d in era['date']])
era['month'] = np.asarray([d.month for d in era['date']])
era['year'] = np.asarray([d.year for d in era['date']])

day = 28

aug_index = np.where(era['day']==day)[0]


msl = nc.variables['msl'][:]
tcwv = nc.variables['tcwv'][:]
uwind = nc.variables['u10'][:]
vwind = nc.variables['v10'][:]
dewpt = nc.variables['d2m'][:]

u = nc2.variables['u'][:]
v = nc2.variables['v'][:]
z = nc2.variables['z'][:]


def wind_magnitude(ucomp,vcomp):
    wind = np.sqrt((ucomp**2) + (vcomp**2))
    return wind


total_wind = wind_magnitude(u,v)
harvey_850wind = total_wind[aug_index,0,:,:]




tcwv_harvey = tcwv[aug_index,:,:]
##Take the daily mean
tcwv_mean = np.nanmean(tcwv_harvey, axis = 0)



#%%

###Plot total column water vapor

plt.figure(figsize=(10,6))
cmin = -12; cmax = 12.; cint = 2.; clevs = np.arange(cmin,cmax,cint)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
xlim = np.array([-99.,-93.]); ylim = np.array([27.,33.])
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
cs = m.contourf(lon2,lat2,tcwv_mean.T,cmap='BuGn',extend='both') 
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('$Kg m^-2$',weight='bold',name='Calibri',size=14)
cticks = []
plt.title('Mean Total Column Water Vapor 8/25/2017',name='Calibri',weight='bold',size=16)
plt.show(block=False)













