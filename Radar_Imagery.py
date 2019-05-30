#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 15:13:03 2019

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np

from scipy import io
import glob
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import sys,getopt
from matplotlib.colors import LinearSegmentedColormap
import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
import pyart



files = glob.glob('*.nc')

sample_file = files[0]
nc = Dataset(sample_file, 'r')

latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
azimuth_data = nc.variables['azimuth'][:]
range_data = nc.variables['gate'][:]


zh = np.ones((360,230,22))*np.nan

for n,j in enumerate(files):
    data = Dataset(j)
    #print(j)
    zh[:,:,n] = data.variables['BaseReflectivity'][:]
    

x = range_data*np.sin(np.deg2rad(azimuth_data))[:,None]
y = range_data*np.cos(np.deg2rad(azimuth_data))[:,None]


#Shape of 360,230,22


    
radar_data_refl = np.ma.array(zh, mask = np.isnan(zh))
proj = cartopy.crs.LambertConformal(central_longitude = nc.RadarLongitude, central_latitude = nc.RadarLatitude)
    
state_borders = cartopy.feature.NaturalEarthFeature(category='cultural', name = 'admin_1_states_provinces_lakes', scale = '50m', facecolor = 'none')

reader = shpreader.Reader('countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())    

#%%

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = 0.; cmax = 80.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)



mesh = ax.pcolormesh(x,y,radar_data_refl[:,:,21], cmap = 'pyart_NWSRef', vmin = 0, vmax = 75)

ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1

ax.set_extent([nc.RadarLongitude - distance_in_degrees, nc.RadarLongitude + distance_in_degrees, nc.RadarLatitude - distance_in_degrees, nc.RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel('dBZ',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
    
  
#ax.plot(-105.1008, 40.1792, 'r',transform=proj)  
    
x = -105.1008
y = 40.1792    
    
geodetic = ccrs.Geodetic()
x,y=proj.transform_point(x,y,geodetic)
#ax.plot(x,y,'ko',markersize=7)  
ax.plot(x,y,'or',markerfacecolor='w')  
   
    
plt.title(r'KFTG $0.5^o$ Base Reflectivity 5/22 17:59 MDT', size = 20) 
plt.savefig('/Users/noahbrauer/Desktop/Consulting/Radar_Images/radar21.png')

plt.show() 
#%%

import imageio

png_dir = '/Users/noahbrauer/Desktop/Consulting/Radar_Images'
images = []
for file_name in os.listdir(png_dir):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('radar.gif', images, duration = 0.5)




