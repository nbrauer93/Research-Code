#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:20:56 2019

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter



import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

file = 'pv_pres.nc'
era = {}

nc = Dataset(file, 'r', unpack = True)

lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:] - 360
lat2, lon2 = np.meshgrid(lat,lon)[:]

level = nc.variables['level'][:]*100
level_mb = level/100

time = nc.variables['time'][:]
timeUnits = nc.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
era['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
era['day'] = np.asarray([d.day for d in era['date']])
era['month'] = np.asarray([d.month for d in era['date']])
era['year'] = np.asarray([d.year for d in era['date']])


pv = nc.variables['pv'][:]
temp = nc.variables['t'][:]



aug_index = np.where((era['month']==8)&(era['day']>24))[0]

pv_aug = pv[aug_index,:,:,:]
temp_aug = temp[aug_index,:,:,:]



#Calculate potential temperature from temperature

def potential_temp(temp,pres):
    theta = temp*(100000/pres)**(287/1004)
    return theta




T,Z,I,J = temp_aug.shape
tempsquish = temp_aug.reshape(Z,T*I*J, order = 'F')

theta_squished = np.ones((28,809760))*np.nan

for i in range(theta_squished.shape[0]):
    for j in range(theta_squished.shape[1]):
        theta_squished[i,j] = potential_temp(tempsquish[i,j],level[i])

#Now reshape back into oringal form
        
theta = theta_squished.reshape(T,Z,I,J)        


#Define a constant latitude and take cross section along this; varies by longitude (-97.5 -94.5 )
#constant latitude is 30degN

theta_hou = theta[:,:,80,:]
pv_hou = pv_aug[:,:,80,:]

theta_hou2 = theta_hou[3,:,350:354]
pv_hou2 = pv_hou[3,:,350:354]

lon_hou = lon[350:354]
#%%


fig = plt.figure(figsize = [10,10])








