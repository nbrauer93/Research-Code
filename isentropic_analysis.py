#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 17:04:00 2018

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

file = 'isentropic.nc'



era = {}

nc = Dataset(file, 'r', unpack = True)



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



aug_index = np.where(era['day']>25)[0]



pres = nc.variables['pres'][:] #350,330,315,300,285
q = nc.variables['q'][:]
mont = nc.variables['mont'][:]
u = nc.variables['u'][:]*1.944
v = nc.variables['v'][:]*1.944

#315 theta surface

thta_315 = pres[aug_index,2,:,:]
q_315 = q[aug_index,2,:,:]*1000

### Wind barbs now at 315K

u_315 = u[aug_index,2,:,:]
v_315 = v[aug_index,2,:,:]

#Do the same at 330K
u_330 = u[aug_index,1,:,:]
v_330 = v[aug_index,1,:,:]

thta_330 = pres[aug_index,1,:,:]
q_330 = q[aug_index,1,:,:]*1000

### Wind barbs now at 330K

u_350 = u[aug_index,1,:,:]
v_350 = v[aug_index,1,:,:]

###Now for  350K

thta_350 = pres[aug_index,0,:,:]
q_350 = q[aug_index,0,:,:]*1000








#%%

plt.figure(figsize=(16,8))
cmin = 0; cmax = 0.00001; cint = 0.000001; clevs = np.arange(cmin,cmax,cint)
#nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
cs = m.contour(lon2,lat2,thta_315[0,:,:].T/100,colors = 'k',extend='both') 
cs2 = m.contourf(lon2,lat2,q_315[0,:,:].T,cmap = 'Greens',extend='both') 
plt.barbs(lon2[::3,::3], lat2[::3,::3], u_315[0,::3,::3], v_315[0,::3,::3], np.sqrt(u_315[0,::3,::3]*u_315[0,::3,::3] + v_315[0,::3,::3]*v_315[0,::3,::3]), fill_empty=True, rounding=False, flagcolor = 'k',
         barbcolor = 'k', sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3))
#m.drawcounties()
plt.clabel(cs, inline=1, fontsize=12, fmt='%d')
CB = plt.colorbar(cs2,shrink=0.5, extend='both', label = '$g kg^-1$')
CB.set_clim(1,15)


#cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1']) 


#cbar = m.colorbar(cs,size='2%')
#cbar.ax.set_ylabel('$m$',weight='bold',name='Calibri',size=14)
#cticks = []
plt.title('300 K Pressure, winds, specific humidity (shaded),  8/26 00Z',name='Calibri',weight='bold',size=16)
plt.show(block=False)




#%%

title = ['8/26 00Z', '8/26 06Z','8/26 12Z','8/26 18Z','8/27 00Z','8/27 06Z','8/27 12Z','8/27 18Z','8/28 00Z','8/28 06Z','8/28 12Z', '8/28 18Z','8/29 00Z', '8/29 06Z', '8/29 12Z','8/29 18Z','8/30 00Z', '8/30 06Z','8/30 12Z', '8/30 18Z', '8/31 00Z', '8/31 06Z', '8/31 12Z', '8/31 18Z']
alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x']
'''
for i in range(u_315.shape[0]):
    plt.figure(figsize=(16,8))
    cmin = 0; cmax = 0.00001; cint = 0.000001; clevs = np.arange(cmin,cmax,cint)
    #nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    cs = m.contour(lon2,lat2,thta_315[i,:,:].T/100,colors = 'k',extend='both') 
    cs2 = m.contourf(lon2,lat2,q_315[i,:,:].T,cmap = 'Greens',extend='both') 
    plt.barbs(lon2[::3,::3], lat2[::3,::3], u_315[i,::3,::3], v_315[i,::3,::3], np.sqrt(u_315[i,::3,::3]*u_315[i,::3,::3] + v_315[i,::3,::3]*v_315[i,::3,::3]), fill_empty=True, rounding=False, flagcolor = 'k',
         barbcolor = 'k', sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3))
    #m.drawcounties()
    plt.clabel(cs, inline=1, fontsize=12, fmt='%d')
    CB = plt.colorbar(cs2,shrink=0.5, extend='both', label = '$g kg^-1$')
    CB.set_clim(1,15)

    plt.title('315 K Pressure, winds, specific humidity (shaded)' + ' ' + str(title[i]) ,name='Calibri',weight='bold',size=16)
    plt.savefig('315K'+str(alphabet[i])+'.png')
    plt.show(block=False)
    '''
#%%    
###330 K surface

for i in range(u_330.shape[0]):
    plt.figure(figsize=(16,8))
    cmin = 0; cmax = 0.00001; cint = 0.000001; clevs = np.arange(cmin,cmax,cint)
    #nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    cs = m.contour(lon2,lat2,thta_330[i,:,:].T/100,colors = 'k',extend='both') 
    cs2 = m.contourf(lon2,lat2,q_330[i,:,:].T,cmap = 'Greens',extend='both') 
    plt.barbs(lon2[::3,::3], lat2[::3,::3], u_330[i,::3,::3], v_330[i,::3,::3], np.sqrt(u_330[i,::3,::3]*u_330[i,::3,::3] + v_330[i,::3,::3]*v_330[i,::3,::3]), fill_empty=False, rounding=False, flagcolor = 'k',
         barbcolor = 'k', sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3))
    #m.drawcounties()
    plt.clabel(cs, inline=1, fontsize=12, fmt='%d')
    CB = plt.colorbar(cs2,shrink=0.5, extend='both', label = '$g kg^-1$')
    CB.set_clim(1,15)

    plt.title('330 K Pressure, winds, specific humidity (shaded)' + ' ' + str(title[i]) ,name='Calibri',weight='bold',size=16)
    #plt.savefig('330k'+str(alphabet[i])+'.png')
    plt.show(block=False)
    



#%%
    
#350 K theta surface


for i in range(u_330.shape[0]):
    plt.figure(figsize=(16,8))
    cmin = 0; cmax = 0.00001; cint = 0.000001; clevs = np.arange(cmin,cmax,cint)
    #nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    xlim = np.array([-120.,-75.]); ylim = np.array([23.,45.])
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    cs = m.contour(lon2,lat2,thta_350[i,:,:].T/100,colors = 'k',extend='both') 
    cs2 = m.contourf(lon2,lat2,q_350[i,:,:].T,cmap = 'Greens',extend='both') 
    plt.barbs(lon2[::3,::3], lat2[::3,::3], u_350[i,::3,::3], v_350[i,::3,::3], np.sqrt(u_350[i,::3,::3]*u_350[i,::3,::3] + v_350[i,::3,::3]*v_350[i,::3,::3]), fill_empty=False, rounding=False, flagcolor = 'k',
         barbcolor = 'k', sizes=dict(emptybarb=0.25, spacing=0.2, height=0.3))
    #m.drawcounties()
    plt.clabel(cs, inline=1, fontsize=12, fmt='%d')
    CB = plt.colorbar(cs2,shrink=0.5, extend='both', label = '$g kg^-1$')
    CB.set_clim(1,15)

    plt.title('350 K Pressure, winds, specific humidity (shaded)' + ' ' + str(title[i]) ,name='Calibri',weight='bold',size=16)
    plt.savefig('350k'+str(alphabet[i])+'.png')
    plt.show(block=False)
        
    
    



#%%
    
###create a  gif here
    
import imageio

png_dir = '/Users/noahbrauer/Desktop/Research/Reanalysis/Isentropic'
images = []
for file_name in os.listdir(png_dir):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('../350KPres.gif', images, duration = 0.5)
    

 


   





    

    
    
    
    
    
    