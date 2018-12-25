#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 16:07:49 2018

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


aug_index = np.where((narr['day']>=24)&(narr['month']==8))[0]


#Extract lhf variable

lhf = files.variables['lhtfl'][aug_index,:,:]

#Calculate the mean and std dev for 1979-2017 period

lhf_mean = np.nanmean(lhf,axis = 0)
lhf_std = np.nanstd(lhf,axis = 0)


lhf_mean_daily25 = np.ones((277,349))*np.nan
lhf_mean_daily26 = np.ones((277,349))*np.nan
lhf_mean_daily27 = np.ones((277,349))*np.nan
lhf_mean_daily28 = np.ones((277,349))*np.nan
lhf_mean_daily29 = np.ones((277,349))*np.nan
lhf_mean_daily30 = np.ones((277,349))*np.nan

for j in range(0,8):
    lhf_mean_daily25[:,:] = np.nanmean(lhf[8*j,:,:], axis = 0)
    lhf_mean_daily26[:,:] = np.nanmean(lhf[(8*j)+1,:,:], axis = 0)
    lhf_mean_daily27[:,:] = np.nanmean(lhf[(8*j)+2,:,:], axis = 0)
    lhf_mean_daily28[:,:] = np.nanmean(lhf[(8*j)+3,:,:], axis = 0)
    lhf_mean_daily29[:,:] = np.nanmean(lhf[(8*j)+4,:,:], axis = 0)
    lhf_mean_daily30[:,:] = np.nanmean(lhf[(8*j)+5,:,:], axis = 0)



#Extract times for Harvey
#%%

file = 'lhtfl.2017.nc'
nc = Dataset(file,'r')
harvey = {}


time2 = nc.variables['time'][:]
timeUnits2 = nc.variables['time'].units
tmpDates2 = num2date(time2,timeUnits2,calendar='gregorian')
harvey['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates2])
harvey['day'] = np.asarray([d.day for d in harvey['date']])
harvey['month'] = np.asarray([d.month for d in harvey['date']])
harvey['year'] = np.asarray([d.year for d in harvey['date']])


harvey_lhf = nc.variables['lhtfl'][:]
harvey_index = np.where((harvey['month']==8)&(harvey['day']>=24))[0]

lhf_harvey = harvey_lhf[harvey_index,:,:]

#Now compute the anomalies and standardize:


lhf_ano = np.ones((8,277,349))

for i in range(lhf_ano.shape[0]):
    lhf_ano[i,:,:] = (lhf_harvey[i,:,:] - lhf_mean)/lhf_std


#Calculate the t statistic for each data point
#Using the 95% confidence interval
#Null hypothesis: The true mean is between lhf_mean +/- the change in LHF (how much it varies); Alternative: Outside of this range

t95 = 1.968  #T-table critical value at 95% confidence interval

lower_bound = np.ones((277,349))*np.nan
upper_bound = np.ones((277,349))*np.nan

for i in range(lower_bound.shape[0]):
    for j in range(lower_bound.shape[1]):
        lower_bound[i,j] = lhf_mean[i,j] - t95*(lhf_std[i,j])/np.sqrt(312-1)
        upper_bound[i,j] = lhf_mean[i,j] + t95*(lhf_std[i,j])/np.sqrt(312-1)
        
    
significant = np.ones((8,277,349))*np.nan

for i in range(significant.shape[0]):
    for j in range(significant.shape[1]):
        for k in range(significant.shape[2]):
            if lhf_harvey[i,j,k]<lower_bound[j,k] or lhf_harvey[i,j,k]>upper_bound[j,k]:
                significant[i,j,k] = lhf_harvey[i,j,k]
            else:
                significant[i,j,k] = np.nan
                





    
    
    
#%%    
    
##Now we plot!


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
cs = m.contourf(lon,lat,lhf_ano[4,:,:],clevs,cmap='RdBu',extend='both') 
m.drawcounties()
m.pcolor(lon, lat, significant[4,:,:], hatch='.', alpha = 0.)
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

plt.title('Latent Heat Flux Anomalies 8/28',name='Calibri',size=14, weight = 'bold')
plt.show(block=False)    




