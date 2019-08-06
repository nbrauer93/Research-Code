#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 15:39:28 2019

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
import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader



import conda
import os
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap





#file = '20170827-033207.netcdf'
file = '20170827-040115.netcdf'
nc = Dataset(file, 'r')


lat = np.linspace(55.0005,20.0005,7000)     
lon = np.linspace(-130.0005,-60.0005,14000)
xlim = [-99,-93]; ylim = [26,32]
ilat = np.where((lat>=ylim[0])&(lat<=ylim[1]))[0]
ilon = np.where((lon>=xlim[0])&(lon<=xlim[1]))[0]
lat_hou = lat[ilat]
lon_hou = lon[ilon]
lathou2,lonhou2 = np.meshgrid(lat_hou,lon_hou)


rotation = nc.variables['RotationTrack60min'][ilat,ilon]


#Mask out values less than 0.01 s^-1


rotation_nan = np.ones((1200,1200))*np.nan

for i in range(rotation.shape[0]):
    for j in range(rotation.shape[1]):
        
        if rotation[i,j]>0.0001:
            rotation_nan[i,j] = rotation[i,j]
            
        else:
            rotation_nan[i,j] = np.nan


#Scale (multiply by 1000)
            
            
rotation_scaled = rotation_nan*1000           







plt.figure(figsize=(10,10))
   
xlim = np.array([-97.5,-94]); ylim = np.array([28.5,31])

cmin = 0.; cmax = 20.; cint = 1.; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='YlOrRd',lut=nlevs)

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lonhou2,lathou2,rotation_scaled.T,clevs,cmap='YlOrRd',extend='both') 

m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel(r'$10^{-3} s^{-1}$',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
    
    
plt.title('3-6 km Azimuthal Shear 0301-0401 UTC 8/27', name = 'Calibri', size = 26)
plt.show(block=False)         






