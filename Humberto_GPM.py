#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:11:25 2019

@author: noahbrauer
"""

import h5py 
import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import pyart



file = '2A.GPM.DPR.V8-20180723.20190915-S160459-E173733.031518.V06A.HDF5'

DPR = h5py.File(file, 'r')

lat2 = DPR['NS']['Latitude'][:,:]    
lon2 = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:]  #nray, nscan
precip_rate = DPR['NS']['SLV']['precipRateNearSurface'][:]
alt = DPR['NS']['PRE']['elevation'][:] #in meters


#set up lat-lon

ind1 = np.where((lon2[:,0]>=-83))
ind2 = np.where((lon2[:,0])<=-70)
ind3 = np.intersect1d(ind1,ind2)



x2 = 2.*17 #48 degrees
re = 6378 #radius of earth in km

theta = -1 *(x2/2.) + (x2/48.)*np.arange(0,49) #Get equal degree increments

theta2 = np.zeros(theta.shape[0]+1)
theta = theta - 0.70833333/2.
theta2[:-1] = theta
theta2[-1] = theta[-1] + 0.70833333

theta = theta2 * (np.pi/180.) #convert to radians

prh = np.zeros([177,50]) #set up matrix 
for i in np.arange(0,177): #loop over num range gates
    for j in np.arange(0,50): #loop over scans 
        a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #407 km is the orbit height 
        prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a) #more geometry 
        
h2 = prh

h3 = np.zeros([h2.shape[0],h2.shape[1]])
for i in np.arange(0,h3.shape[1]):
    h3[:,i] = h2[::-1,i] #reverse order so as index go up the height goes up
    
    
from pyproj import Proj

#Extract the vertical profile of reflectivity

ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:]  #176 levels

ku = ku[:,27,:]
lons = DPR['NS']['Longitude'][ind3,27]
lats = DPR['NS']['Latitude'][ind3,27]

# pick the starting point to be the reference point to calc dist. 
lat0 = lats[0]
lon0 = lons[0]
#this projection is not exactly right, but works 
p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0)


lat_3d = np.zeros(ku.shape)
lon_3d = np.zeros(ku.shape)
for i in np.arange(0,ku.shape[0]):
    lat_3d[i,:] = lats[i]
    lon_3d[i,:] = lons[i]
#convert to distances 
x,y = p(lon_3d,lat_3d)

R_gpm = np.sqrt(x**2  + y**2)*np.sign(x)

ku = ku[:,::-1]


#remove bad data
ku = np.ma.masked_where(ku <= 12,ku)


y = np.zeros([ku.shape[0],ku.shape[1]])


h4 = h3[:,27]
    
for i in np.arange(y.shape[1]):
    y[:,i] = h4[i]


plt.figure(figsize=(10,10))

vmax = 60
vmin = 12

R_min = R_gpm.min()
R_max = R_gpm.max()

pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
plt.xlabel('N-S Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title('GPM Overpass 9/15 1604 UTC N-S Cross-Section', size = 20)
plt.xlim(600,1500)
plt.ylim(0,15)
plt.show()


















'''


#Filter out zeros and convert to NaN

z_nan = np.ones((7936,49))*np.nan
precip_nan = np.ones((7936,49))*np.nan

for i in range(z.shape[0]):
    for j in range(z.shape[1]):
        
        if z[i,j] <0 or precip_rate[i,j] >100:
            z_nan[i,j] = np.nan
            precip_nan[i,j] = np.nan
            
        else:
            z_nan[i,j] = z[i,j]
            precip_nan[i,j] = precip_rate[i,j]
            
            
            





#Plot radar reflectivity


cmin = 0.; cmax = 75.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  

#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-70]); ylim = np.array([26,35])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lon2,lat2,z_nan,clevs,cmap='pyart_NWSRef',extend='both') 


#m.drawcounties()

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dBZ',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
   
plt.title('GPM Overpass 9/15/2019 1604 UTC', size = 20)




plt.show(block=False)



#plot precipitation rate

cmin = 0.; cmax = 100.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  

#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-70]); ylim = np.array([26,35])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lon2,lat2,precip_nan,clevs,cmap='pyart_NWSRef',extend='both') 


#m.drawcounties()

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('mm/hr',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
   
plt.title('GPM Precipitation Rate 9/15/2019 1604 UTC', size = 20)




plt.show(block=False)

'''

