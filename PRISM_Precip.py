#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 11:43:56 2018

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from mpl_toolkits.basemap import Basemap,maskoceans,interp,shiftgrid
from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob




#files = 'PRISM_ppt_stable_4kmD2_20170828_bil.bil'

prism_path = glob.glob('*_bil.bil')

 
def read_prism_hdr(hdr_path):
    """Read an ESRI BIL HDR file"""
    with open(hdr_path, 'r') as input_f:
        header_list = input_f.readlines()
    return dict(item.strip().split() for item in header_list)
 
def read_prism_bil(bil_path):
    """Read an array from ESRI BIL raster file"""
    hdr_dict = read_prism_hdr(bil_path.replace('.bil', '.hdr'))
    # For now, only use NROWS, NCOLS, and NODATA
    # Eventually use NBANDS, BYTEORDER, LAYOUT, PIXELTYPE, NBITS
 
    prism_array = np.fromfile(bil_path, dtype=np.float32)
    prism_array = prism_array.reshape(
        int(hdr_dict['NROWS']), int(hdr_dict['NCOLS']))
    prism_array[prism_array == float(hdr_dict['NODATA'])] = np.nan
    return prism_array
 
#prism_path = 'PRISM_ppt_stable_4kmD2_20170828_bil.bil'
    
prism_array = np.ones((621,1405,5))*np.nan

for i in range(len(prism_path)):
    prism_array[:,:,i] = read_prism_bil(prism_path[i])

#prism_array = read_prism_bil(prism_path)


# making the data grid points for PRISM

xlim = np.array([-125.0208333, -66.4791667])
ylim = np.array([24.0625000, 49.9375000])
prism_cols = 1405
prism_rows = 621

x = np.linspace(-125.0208333, -66.4791667, 1405)
y = np.linspace(24.0625000, 49.9375000, 621)[::-1]
 
iX = np.where( (x>=xlim[0]) & (x<=xlim[1]) )[0]
iY = np.where( (y>=ylim[0]) & (y<=ylim[1]) )[0]
prism_nodata = -9999
lat,lon = np.meshgrid(y[iY],x[iX])


total_precip  = np.nansum(prism_array, axis = 2)




#%%

date = ['8/26','8/27','8/28','8/29','8/30']


for i in range(prism_array.shape[0]-1):
    plt.figure(figsize=(10,6))
    cmin = 0.; cmax = 500.; cint = 25; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-97,-93]); ylim = np.array([28,31])

    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates()

    cs = m.contourf(lon,lat,prism_array[:,:,i].T,clevs,cmap='GnBu',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
    m.drawcounties()
    cbar = m.colorbar(cs,size='2%')
    cbar.ax.set_ylabel('mm',weight='bold',name='Calibri',size=14)
    cticks = []
    for i in clevs:
        cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
    for i in cbar.ax.yaxis.get_ticklabels():
        i.set_family('Calibri')
        i.set_size(14)
    plt.title('Precipitation',name='Calibri',weight='bold',size=16)
    plt.show(block="False")




#%%



plt.figure(figsize=(10,6))
cmin = 0.; cmax = 1400.; cint = 50; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
xlim = np.array([-97,-93]); ylim = np.array([28,31])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='h')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates()

cs = m.contourf(lon,lat,total_precip.T,clevs,cmap='GnBu',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('mm',weight='bold',name='Calibri',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
    
    
x2star,y2star = m(-95.3698,29.7604)
m.plot(x2star,y2star,'ro',markersize=7)
label = 'Houston'
plt.text(x2star+0.1,y2star+0.1,label)
        
plt.title('Total Precipitation 8/26-8/30 2017',name='Calibri',weight='bold',size=16)
plt.show(block="False")


