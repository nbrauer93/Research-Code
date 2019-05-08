#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:03:47 2019

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

file = 'KHGX_SDUS54_N0VHGX_201708270331.nc'
file2 = 'KHGX_SDUS54_N0RHGX_201708270331.nc'
file3 = 'KHGX_SDUS24_N1UHGX_201708270331.nc'
file4 = 'KHGX_SDUS24_N1SHGX_201708270331.nc'
file5 = 'KHGX_SDUS54_N0SHGX_201708270331.nc'


nc = Dataset(file, 'r')
nc2 = Dataset(file2, 'r')
nc3 = Dataset(file3, 'r')
nc4 = Dataset(file4, 'r')
nc5 = Dataset(file5, 'r')

latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
velocity = nc.variables['RadialVelocity'][:]
range_data = nc.variables['gate'][:]
azimuth_data = nc.variables['azimuth'][:]
#0.5deg reflectivity
zh = nc2.variables['BaseReflectivity'][:]
#1.5deg velocity
velocity15 = nc3.variables['BaseVelocityDV'][:]
azimuth_data15 = nc3.variables['azimuth'][:]
range_data15 = nc3.variables['gate'][:]

#1.5deg storm-relative
vel_sr = nc4.variables['StormMeanVelocity'][:]
azimuth_datasr = nc4.variables['azimuth'][:]
range_datasr = nc4.variables['gate'][:]


#0.5 deg Storm-relative
vel_sr2 = nc5.variables['StormMeanVelocity'][:]
azimuth_datasr2 = nc5.variables['azimuth'][:]
range_datasr2 = nc5.variables['gate'][:]



x = range_data*np.sin(np.deg2rad(azimuth_data))[:,None]
y = range_data*np.cos(np.deg2rad(azimuth_data))[:,None]

x2 = range_data15*np.sin(np.deg2rad(azimuth_data15))[:,None]
y2 = range_data15*np.cos(np.deg2rad(azimuth_data15))[:,None]

x3 = range_datasr*np.sin(np.deg2rad(azimuth_datasr))[:,None]
y3 = range_datasr*np.cos(np.deg2rad(azimuth_datasr))[:,None]

x3 = range_datasr*np.sin(np.deg2rad(azimuth_datasr))[:,None]
y3 = range_datasr*np.cos(np.deg2rad(azimuth_datasr))[:,None]

x4 = range_datasr2*np.sin(np.deg2rad(azimuth_datasr2))[:,None]
y4 = range_datasr2*np.cos(np.deg2rad(azimuth_datasr2))[:,None]


radar_data = np.ma.array(velocity, mask = np.isnan(velocity))
radar_data_refl = np.ma.array(zh, mask = np.isnan(zh))
radar_data_vel = np.ma.array(velocity15, mask = np.isnan(velocity15))
radar_data_srvel = np.ma.array(vel_sr, mask = np.isnan(vel_sr)) 
radar_data_srvel2 = np.ma.array(vel_sr2, mask = np.isnan(vel_sr2)) 


proj = cartopy.crs.LambertConformal(central_longitude = nc.RadarLongitude, central_latitude = nc.RadarLatitude)
proj15 = cartopy.crs.LambertConformal(central_longitude = nc3.RadarLongitude, central_latitude = nc3.RadarLatitude)
projsr = cartopy.crs.LambertConformal(central_longitude = nc4.RadarLongitude, central_latitude = nc4.RadarLatitude)
projsr2 = cartopy.crs.LambertConformal(central_longitude = nc5.RadarLongitude, central_latitude = nc5.RadarLatitude)




state_borders = cartopy.feature.NaturalEarthFeature(category='cultural', name = 'admin_1_states_provinces_lakes', scale = '50m', facecolor = 'none')

reader = shpreader.Reader('countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = -75.; cmax = 75.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSVel',lut=nlevs)



mesh = ax.pcolormesh(x,y,radar_data, cmap = 'pyart_NWSVel' )
ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1.

ax.set_extent([nc.RadarLongitude - distance_in_degrees, nc.RadarLongitude + distance_in_degrees, nc.RadarLatitude - distance_in_degrees, nc.RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel(r'$m s^{-1} $',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('KHGX Radial Velocity 8/27 0331 UTC', size = 24)    


#1.5 deg tilt

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = -75.; cmax = 75.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSVel',lut=nlevs)



mesh = ax.pcolormesh(x2,y2,radar_data_vel, cmap = 'pyart_NWSVel' )
ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1.

ax.set_extent([nc3.RadarLongitude - distance_in_degrees, nc3.RadarLongitude + distance_in_degrees, nc3.RadarLatitude - distance_in_degrees, nc3.RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel(r'$m s^{-1} $',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('KHGX 1.5 deg Radial Velocity  8/27 0331 UTC', size = 24)  



#1.5deg S-R velocity


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = -75.; cmax = 75.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSVel',lut=nlevs)



mesh = ax.pcolormesh(x3,y3,radar_data_srvel, cmap = 'pyart_NWSVel', vmin = -45, vmax = 45  )
ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1.

ax.set_extent([nc3.RadarLongitude - distance_in_degrees, nc3.RadarLongitude + distance_in_degrees, nc3.RadarLatitude - distance_in_degrees, nc3.RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel(r'$m s^{-1} $',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('KHGX 1.5 deg Storm-Relative Velocity  8/27 0331 UTC', size = 24)  

###0.5deg S-R velocity



fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = -75.; cmax = 75.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSVel',lut=nlevs)



mesh = ax.pcolormesh(x4,y4,radar_data_srvel2, cmap = 'pyart_NWSVel', vmin = -45, vmax = 45 )
ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1.

ax.set_extent([nc3.RadarLongitude - distance_in_degrees, nc3.RadarLongitude + distance_in_degrees, nc3.RadarLatitude - distance_in_degrees, nc3.RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel(r'$m s^{-1} $',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('KHGX 0.5 deg Storm-Relative Velocity  8/27 0331 UTC', size = 24)  




#####Plot base reflectivity

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = 0.; cmax = 80.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)



mesh = ax.pcolormesh(x,y,radar_data_refl, cmap = 'pyart_NWSRef', vmin = 0, vmax = 75 )
ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray')

distance_in_degrees = 1.

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
plt.title('KHGX Base Reflectivity 8/27 0331 UTC', size = 24)    




