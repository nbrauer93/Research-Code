#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 11:50:03 2019

@author: noahbrauer
"""

#Script for basic plotting of GPM HDF5 data

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

import h5py 
import numpy as np
import matplotlib.pyplot as plt

file = '3B-HHR-L.MS.MRG.3IMERG.20180914-S010000-E012959.0060.V05B.HDF5'


f = h5py.File(file, 'r')
for key in f.keys():
    print(key) #MS HS
    
group = f[key]   


#Print keys inside group
for key in group.keys():
    print(key)
    
#Extract data

lat = group['lon'].value     
lon = group['lat'].value
precip = group['HQprecipitation'].value


#Define a lat-lon grid


lat2,lon2 = np.meshgrid(lat,lon)[:]


#Make nans


precip_nan = np.ones((3600,1800))*np.nan

for i in range(precip.shape[0]):
    for j in range(precip.shape[1]):
        
        if precip[i,j] == 0:
            precip_nan[i,j] = np.nan
        else:
            precip_nan[i,j] = precip[i,j]


#Plot

#%%

cmin = 0; cmax = 20; cint = 0.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='GnBu',lut=nlevs)
   
plt.figure(figsize=(10,6))

xlim = np.array([-30,20]); ylim = np.array([40,70])

m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m.contourf(lon2,lat2,precip_nan.T,clevs,cmap='GnBu',extend='both') 

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('mm',weight='bold',size=16)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)

plt.title('GPM 1.5 Hour Precipitation',name='Calibri',size=16, weight = 'bold')
plt.show(block=False) 



#%%


file_dpr = 'GPMCOR_KAR_1908261440_1612_031206_1BS_DAB_05C.h5'

f2 = h5py.File(file_dpr, 'r')


for key in f2.keys():
    print(key)
    
group_dpr = f2[key]    

for key in group_dpr.keys():
    print(key)
        
    
vert = np.array(f2['MS']['VertLocate']) 

#Read in attributes






    

 

    










    
    
    
    
    
    




 

