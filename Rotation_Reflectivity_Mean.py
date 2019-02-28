#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 15:20:20 2019

@author: noahbrauer
"""

from Temporal_1hr_mean import read_file, diff_reflect,reflect_ncdc
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob
from matplotlib.colors import LinearSegmentedColormap
import warnings
import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap



files_rot = glob.glob('*.netcdf')
files_radar = glob.glob('*.nc')


radar_file = 'nexrad_3d_v4_0_20170827T222000Z.nc'
nc = Dataset(radar_file,'r')
radarlat = nc.variables['Latitude'][:]
radarlon = nc.variables['Longitude'][:] -360
lat2,lon2 = np.meshgrid(radarlat,radarlon)


zh = np.ones((288,528,672))*np.nan
zdr = np.ones((288,528,672))*np.nan
kdp = np.ones((288,528,672))*np.nan


for n,j in enumerate(files_radar):
    data = read_file(j)
    print(j)
    zdr[n,:,:] = data['zdr']['values'][4,:,:]
    kdp[n,:,:] = data['kdp']['values'][4,:,:]
    zh[n,:,:] = data['Z_H']['values'][4,:,:]


#Take 1-hour means of reflectivity, zdr, and  kdp


zhmean = np.ones((24,528,672))*np.nan
zdrmean = np.ones((24,528,672))*np.nan
kdpmean = np.ones((24,528,672))*np.nan

for i in range(zhmean.shape[0]):
    zhmean[i,:,:] = np.nanmean(zh[12*i:(12*i)+12,:,:], axis = 0)
    zdrmean[i,:,:] = np.nanmean(zdr[12*i:(12*i)+12,:,:], axis = 0)
    kdpmean[i,:,:] = np.nanmean(kdp[12*i:(12*i)+12,:,:], axis = 0)




lat = np.linspace(55.0005,20.0005,7000)     
lon = np.linspace(-130.0005,-60.0005,14000)

testlat,testlon = np.meshgrid(lat,lon)

xlim = [-99,-93];ylim = [26,32]
ilat = np.where((lat>=ylim[0])&(lat<=ylim[1]))[0]
ilon = np.where((lon>=xlim[0])&(lon<=xlim[1]))[0]
latshrink = lat[ilat]
lonshrink = lon[ilon]
lathou,lonhou = np.meshgrid(latshrink,lonshrink)


rotation = np.ones((1200,1200,26))*np.nan
    
for t,s in enumerate(files_rot):
    rot = Dataset(s)
    print(s)
    rotation[:,:,t] = rot.variables['RotationTrack60min'][ilat,ilon]
    
sig_rotation = np.ones((1200,1200,26))*np.nan
for i in range(rotation.shape[0]):
    for j in range(rotation.shape[1]):
        for k in range(rotation.shape[2]):
            
            if rotation[i,j,k]>0.001:
                sig_rotation[i,j,k] = rotation[i,j,k]


#%%



#add in the colormaps 


zhcolor = reflect_ncdc()
zdrcolor = diff_reflect()

cmin = 0.; cmax = 75.; cint = 5.; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=zhcolor,lut=nlevs)
    
plt.figure(figsize=(10,6))
   
xlim = np.array([-97,-94]); ylim = np.array([29,30.5])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lon2,lat2,zhmean[3,:,:].T,clevs,cmap=zhcolor,extend='both') 
cs2 = m.contour(lonhou.T,lathou.T,sig_rotation[:,:,3], colors = 'k', linewidths = 1.5) 


m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dBZ',weight='bold',name='Calibri',size=14)
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
plt.text(x2star+0.05,y2star+0.05,label)    
      
plt.title('2 km Mean Reflectivity Factor',name='Calibri',weight='bold',size=16)
plt.show(block=False)         







