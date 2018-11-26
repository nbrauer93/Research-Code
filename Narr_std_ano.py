#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 13:07:01 2018

@author: noahbrauer
"""

import matplotlib.pyplot as plt 
from netCDF4 import Dataset, num2date, MFDataset
import netCDF4 as nc
import numpy as np
from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


files = nc.MFDataset('*.nc')

latitude = files.variables['lat'][:]
longitude = files.variables['lon'][:]
pwat = files.variables['pr_wtr'][:]
#%%

xlim = np.array([-100,-85])
ylim = np.array([24,34])

ilon = np.where((longitude>=xlim[0])&(longitude<=xlim[1]))[0]
ilat = np.where((latitude>=ylim[0])&(latitude<=ylim[1]))[0]

narr = {}

time = files.variables['time'][:]
timeUnits = files.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])


aug_index_std = np.where((narr['month']==8))[0]

pwat_aug = pwat[aug_index_std,:,:]

#Standard deviation for pwat for all days in August (one value for each point in space)

aug_std_pwat = np.nanstd(pwat_aug,axis = 0)

#Take the mean for EACH day

day_index = np.where((narr['month']==8)&(narr['day']>=25))[0]
day_pwat = pwat[day_index,:,:]



harvey_mean = np.ones((7,277,349))*np.nan

for i in range(harvey_mean.shape[0]):
    harvey_mean[i,:,:] = np.nanmean(day_pwat,axis=0)


##Now the actual days we are concerned about (7 of them 8/25-8/31)
harvey_index = np.where((narr['month']==8)&(narr['day']>=25)&(narr['year']==2017))[0]
pwat_harvey = pwat[harvey_index,:,:]



#Now standardize your anomalies:


pwat_ano_std = (pwat_harvey-harvey_mean)/aug_std_pwat



#%%

#Now we plot

title_date = ['8/25','8/26','8/27','8/28','8/29','8/30','8/31']

for i in range(pwat_ano_std.shape[0]):

    plt.figure(figsize=(10,6))
    cmin = -4.; cmax = 4.; cint = 0.5; clevs = np.arange(cmin,cmax,cint)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-101.,-93.]); ylim = np.array([25.,33.])
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  
    cs = m.contourf(longitude,latitude,pwat_ano_std[i,:,:],cmap='BrBG',extend='both') 
    m.drawcounties()
    cbar = m.colorbar(cs,size='2%')
    cbar.ax.set_ylabel('$\sigma$',weight='bold',name='Calibri',size=14)
    cticks = []
    
    
    x2star,y2star = m(-95.3698,29.7604)
    m.plot(x2star,y2star,'ro',markersize=7)
    label = 'Houston'
    plt.text(x2star+0.1,y2star+0.1,label)
    
    x3star,y3star = m(-97.3964,27.8006)
    m.plot(x3star,y3star,'ro',markersize=7)
    label2 = 'Corpus Christi'
    plt.text(x3star+0.25,y3star-0.1,label2)

    plt.title('Precipitable Water Anomalies (Standardized)' + ' '  + str(title_date[i]),name='Calibri',weight='bold',size=16, y = 1.05)
    title_string = 'From 1979-2017 Climatology'
    plt.suptitle(title_string, y=0.92, fontsize=12)
    plt.show(block=False) 











