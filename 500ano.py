#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:49:07 2019

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


file = 'geopotential.nc'
filestd = 'interim_stdev_201708.nc'
filemean = 'interim_mean_201708.nc'


nc = Dataset(file, 'r', unpack = True)
nc2 = Dataset(filestd, 'r', unpack = True)
nc3 = Dataset(filemean, 'r', unpack = True)


era = {}
era2 = {}



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



###For the mean and std files:
time2 = nc2.variables['time'][:]
timeUnits2 = nc2.variables['time'].units
tmpDates2 = num2date(time2,timeUnits2,calendar='gregorian')
era2['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates2])
era2['day'] = np.asarray([d.day for d in era2['date']])
era2['month'] = np.asarray([d.month for d in era2['date']])
era2['year'] = np.asarray([d.year for d in era2['date']])

aug_index = np.where(era['day']>21)[0]
aug_index_climo = np.where(era2['day']>21)[0]

zstd = nc2.variables['z'][aug_index_climo,2,:,:]
zmean = nc3.variables['z'][aug_index_climo,2,:,:]
z = nc.variables['z'][aug_index,:,:]



z_ano  = (z - zmean)-zstd


#Now plot it
#%%
title = ['8/22 00Z','8/22 12Z','8/23 00Z','8/23 12Z','8/24 00Z', '8/24 12Z', '8/25 00Z', '8/25 12Z',  '8/26 00Z','8/26 12Z','8/27 00Z','8/27 12Z','8/28 00Z','8/28 12Z','8/29 00Z', '8/29 12Z','8/30 00Z','8/30 12Z', '8/31 00Z', '8/31 12Z']
alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t']

for i in range(z.shape[0]):
    cmin = -6.; cmax = 6.; cint = 1; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    
    plt.figure(figsize=(10,6))
   
    xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])

    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
    m.drawcoastlines(); m.drawstates(), m.drawcountries()  
    cs = m.contourf(lon2,lat2,z_ano[i,:,:].T,clevs,cmap='bwr',extend='both') 
    cs2 = m.contour(lon2.T,lat2.T,z[i,:,:]/100, colors = 'k',extend = 'both') 
    plt.clabel(cs2, inline=1, fontsize=12, fmt='%d')


    
    cbar = m.colorbar(cs,size='2%')
    cbar.ax.set_ylabel(r'$\sigma$',name='Calibri',size=14)
    cticks = []
    for i in clevs:
        cticks.append(int(i)) if i.is_integer() else cticks.append(i)
        cbar.set_ticks(clevs[::4])
        cbar.set_ticklabels(cticks[::4])
    for i in cbar.ax.yaxis.get_ticklabels():
        i.set_family('Calibri')
        i.set_size(14)
        
    x2star,y2star = m(-95.3698,29.7604)
    m.plot(x2star,y2star,'ro',markersize=7, color = 'k')
    label = 'Houston'
    plt.text(x2star+0.05,y2star+0.05,label, weight = 'bold',size = 10)    

    #plt.title('500 mb Geopotential Height (dam) and Standardized Anomalies' + ' '  + str(title[i]),name='Calibri',weight='bold',size=22)
    #plt.savefig('500ano'+str(alphabet[i])+'.png')      

    plt.show(block=False) 





