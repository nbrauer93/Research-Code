#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 14:21:19 2019

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob
from scipy.interpolate import griddata


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

files = MFDataset('*.nc','r') #1979-2017 climatology


narr = {}
lat = files.variables['lat'][:]
lon = files.variables['lon'][:] 

time = files.variables['time'][:]
timeUnits = files.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])


aug_index = np.where((narr['day']>=25)&(narr['month']==8))[0]
harvey_index = np.where((narr['day']>=25)&(narr['month']==8)&(narr['year']==2017))[0]

mconv = files.variables['mconv'][aug_index,:,:]
mconv_harvey = files.variables['mconv'][harvey_index,:,:] #Multiply by -1 to get units of moisture flux convergence


uwind = 'uwnd.10m.2017.nc'
vwind = 'vwnd.10m.2017.nc'


ncu = Dataset(uwind, 'r')
ncv = Dataset(vwind,'r')


u = ncu.variables['uwnd'][236:241,:,:]
v = ncv.variables['vwnd'][236:241,:,:]



files.close()
#%%
mconv = mconv*-1
mconv_harvey = mconv_harvey*-1


#Take each day from Aug 25-31 for all years 1979-2017


mconv25 = np.ones((38,277,349))*np.nan
mconv26 = np.ones((38,277,349))*np.nan
mconv27 = np.ones((38,277,349))*np.nan
mconv28 = np.ones((38,277,349))*np.nan
mconv29 = np.ones((38,277,349))*np.nan
mconv30 = np.ones((38,277,349))*np.nan


for i in range(mconv25.shape[0]):
    mconv25[i,:,:] = mconv[i*7,:,:]
    mconv26[i,:,:] = mconv[(i*7)+1,:,:]
    mconv27[i,:,:] = mconv[(i*7)+2,:,:]
    mconv28[i,:,:] = mconv[(i*7)+3,:,:]
    mconv29[i,:,:] = mconv[(i*7)+4,:,:]
    mconv30[i,:,:] = mconv[(i*7)+5,:,:]
    

###Now take daily means/stds:
    
mean25 = np.nanmean(mconv25, axis = 0)
mean26 = np.nanmean(mconv26, axis = 0)
mean27 = np.nanmean(mconv27, axis = 0)
mean28 = np.nanmean(mconv28, axis = 0)
mean29 = np.nanmean(mconv29, axis = 0)
mean30 = np.nanmean(mconv30, axis = 0)  


std25 = np.nanstd(mconv25, axis = 0)  
std26 = np.nanstd(mconv26, axis = 0) 
std27 = np.nanstd(mconv27, axis = 0) 
std28 = np.nanstd(mconv28, axis = 0) 
std29 = np.nanstd(mconv29, axis = 0) 
std30 = np.nanstd(mconv30, axis = 0) 



#Now standardize the data (daily anomalies)



mconv25_ano = (mconv_harvey[0,:,:] - mean25)/std25 
mconv26_ano = (mconv_harvey[1,:,:] - mean26)/std26 
mconv27_ano = (mconv_harvey[2,:,:] - mean27)/std27 
mconv28_ano = (mconv_harvey[3,:,:] - mean28)/std28 
mconv29_ano = (mconv_harvey[4,:,:] - mean29)/std29 
mconv30_ano = (mconv_harvey[5,:,:] - mean30)/std30 




#%%


cmin = -4.1; cmax = 4.1; cint = 0.2; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,6))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-99,-92]); ylim = np.array([25,32])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,mconv25_ano,clevs,cmap='BrBG',extend='both') 
m.drawcounties()
m.quiver(lon,lat,u[0,:,:], v[0,:,:],scale = 400)
#m.pcolor(lon, lat, significant27, hatch='.', alpha = 0.)
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel(r'[$\sigma$]',weight='bold',size=16)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)

plt.title('Moisture Flux Convergence Anomalies 8/25',name='Calibri',size=16, weight = 'bold')
plt.show(block=False)   






