#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:22:47 2019

@author: noahbrauer
"""


import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np

from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator

import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader


file = '0.5_deg_vel201708270331.nc'
nc = Dataset(file, 'r')

velocity = nc.variables['BaseVelocityDV'][:]
range_data = nc.variables['gate'][:]
azimuth_data = nc.variables['azimuth'][:]

x = range_data*np.sin(np.deg2rad(azimuth_data))[:,None]
y = range_data*np.cos(np.deg2rad(azimuth_data))[:,None]


#radar_data = np.ma.array(velocity, mask = np.isnan(velocity))


u_r = np.diff(velocity,axis = 1)
u_s = np.diff(velocity, axis = 0)


def magnitude(x,y):
    mag = np.sqrt((x**2)+(y**2))
    return mag


az_shear = magnitude(u_r[:-1,:],u_s[:,:-1])

az_data = np.ma.array(az_shear, mask = np.isnan(az_shear))



#Filter out all values less than 5 m/s


'''
az_data_filter = np.copy(az_data)
az_data_filter[az_data_filter >= 5] = np.nan
'''


az_data_filter = np.ones((359,1199))*np.nan

for i in range(az_data.shape[0]):
    for j in range(az_data.shape[1]):
        
        if az_data[i,j]<5:
            az_data_filter[i,j] = np.nan
        else:
            az_data_filter[i,j] = az_data[i,j]



proj = cartopy.crs.LambertConformal(central_longitude = nc.RadarLongitude, central_latitude = nc.RadarLatitude)


state_borders = cartopy.feature.NaturalEarthFeature(category='cultural', name = 'admin_1_states_provinces_lakes', scale = '50m', facecolor = 'none')

reader = shpreader.Reader('countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())




fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = 0.; cmax = 20.; cint = 2; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='Reds',lut=nlevs)



mesh = ax.pcolormesh(x,y,az_data, cmap = 'Reds', vmin = 0, vmax = 15  )
ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1.

ax.set_extent([nc.RadarLongitude - distance_in_degrees, nc.RadarLongitude + distance_in_degrees, nc.RadarLatitude - distance_in_degrees, nc.RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel(r'$ms^{-1}$',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    
    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])

for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)

plt.title(r'KHGX $0.5^o$ Azimuthal Shear 8/27 0331 UTC', size = 20)
plt.show()