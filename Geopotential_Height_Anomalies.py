#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 16:24:29 2018

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

file = 'era.nc'
file2 = 'shear.nc'
file3 = 'interim_mean_201708.nc'
file4 = 'interim_stdev_201708.nc'


era = {}
era2 = {}

nc = Dataset(file, 'r', unpack = True)
nc2 = Dataset(file2, 'r', unpack = True)
nc3 = Dataset(file3, 'r', unpack = True)
nc4 = Dataset(file4, 'r', unpack = True)


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
time2 = nc4.variables['time'][:]
timeUnits2 = nc4.variables['time'].units
tmpDates2 = num2date(time2,timeUnits2,calendar='gregorian')
era2['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates2])
era2['day'] = np.asarray([d.day for d in era2['date']])
era2['month'] = np.asarray([d.month for d in era2['date']])
era2['year'] = np.asarray([d.year for d in era2['date']])

#%%

aug_index = np.where(era['day']>21)[0]
aug_index_climo = np.where(era2['day']>21)[0]



msl = nc.variables['msl'][:]
tcwv = nc.variables['tcwv'][:]
uwind = nc.variables['u10'][:]
vwind = nc.variables['v10'][:]
dewpt = nc.variables['d2m'][:]

u = nc2.variables['u'][:]*1.944   #convert from m/s to knots
v = nc2.variables['v'][:]*1.944
z = nc2.variables['z'][:]

zmean = nc3.variables['z'][:]
zstd = nc4.variables['z'][:]

zmean500 = zmean[:,2,:,:]
zstd500 = zstd[:,2,:,:]

zmean_aug = zmean500[aug_index_climo,:,:]
zstd_aug = zstd500[aug_index_climo,:,:]

z_aug = z[aug_index,1,:,:][::2]
z_ano = (z_aug - zmean_aug)/zstd_aug



#####Plot 250 hpa geopotential height, isotachs

z250 = z[:,0,:,:]
u250 = u[:,0,:,:]
v250 = v[:,0,:,:]

#Take the magnitude of the wind

def wind_magnitude(ucomp,vcomp):
    wind = np.sqrt((ucomp**2) + (vcomp**2))
    return wind


wind250 = wind_magnitude(u250,v250)


wind250_no_nan = wind250[np.logical_not(np.isnan(wind250))]
z250_no_nan = z250[np.logical_not(np.isnan(z250))]
T,I,J = wind250.shape
wind250_nonan = wind250_no_nan.reshape(T,I,J, order = False)
z250_nonan = z250_no_nan.reshape(T,I,J, order = False)





#%%


###Plot here


title = ['8/22 00Z','8/22 12Z','8/23 00Z','8/23 12Z','8/24 00Z', '8/24 12Z', '8/25 00Z', '8/25 12Z',  '8/26 00Z','8/26 12Z','8/27 00Z','8/27 12Z','8/28 00Z','8/28 12Z','8/29 00Z', '8/29 12Z','8/30 00Z','8/30 12Z', '8/31 00Z', '8/31 12Z']
alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t']


for i in range(z_aug.shape[0]):

    plt.figure(figsize=(16,16))
    cmin = 0; cmax = 70; cint = 5.; clevs = np.arange(cmin,cmax,cint)
    #nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    cs = m.contour(lon2,lat2,z_aug[i,:,:].T/100,colors = 'k',extend='both') 
    cs2 = m.contourf(lon2,lat2,z_ano[i,:,:].T,cmap = 'coolwarm',extend='both') 
    #m.drawcounties()
    plt.clabel(cs, inline=1, fontsize=12, fmt='%d')
    CB = plt.colorbar(cs2, shrink=0.4, extend='both')
    CB.ax.set_ylabel('$\sigma$',weight='bold',name='Calibri',size=14)
    CB.set_clim(-2.0, 2.0)

    plt.title('500 mb Geopotential Height (dam) and Standardized Anomalies' + ' '  + str(title[i]),name='Calibri',weight='bold',size=22)
    plt.savefig('500ano'+str(alphabet[i])+'.png')
    
    plt.show(block=False)
    
    #import imageio

'''
png_dir = '/Users/noahbrauer/Desktop/Research/Reanalysis/Heights'
images = []
for file_name in os.listdir(png_dir):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('500heightano.gif', images, duration = 0.5)
'''    


#%%

#Do the same for 250 hPa


for i in range(u250.shape[0]):
    plt.figure(figsize=(10,6))
    cmin = 0; cmax = 0.00001; cint = 0.000001; clevs = np.arange(cmin,cmax,cint)
    #nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    cs = m.contour(lon2[::4,::4],lat2[::4,::4],z250_nonan[i,::4,::4].T/100,colors = 'k',extend='both') 
    cs2 = m.contourf(lon2[::4,::4],lat2[::4,::4],wind250_nonan[i,::4,::4].T,cmap = 'Greens',extend='both') 
    plt.barbs(lon2[::4,::4], lat2[::4,::4], u250[i,::4,::4], v250[i,::4,::4], wind250[i,::4,::4], fill_empty=True, rounding=False, flagcolor = 'k',
         barbcolor = 'k', sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3))
    #m.drawcounties()
    #plt.clabel(cs, inline=1, fontsize=12, fmt='%d')
    CB = plt.colorbar(cs2,shrink=0.5, extend='both', label = 'knots')
    CB.set_clim(0,120)

    plt.title('250 hPa Geopotential Height (dam), Isotachs' + ' ' + str(title[i]) ,name='Calibri',weight='bold',size=16)
    #plt.savefig('330k'+str(alphabet[i])+'.png')
    plt.show(block=False)
    






