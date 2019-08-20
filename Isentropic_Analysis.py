#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 13:58:19 2019

@author: noahbrauer
"""

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import xarray as xr
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, add_timestamp
from metpy.units import units


#List the files


hgt_file = 'hgt.201902.nc'
temp_file = 'air.201902.nc'
q_file = 'shum.201902.nc'
u_file = 'uwnd.201902.nc'
v_file = 'vwnd.201902.nc'

#Read in each file 

nc_hgt = Dataset(hgt_file, 'r')
nc_temp = Dataset(temp_file, 'r')
nc_q = Dataset(q_file, 'r')
nc_u = Dataset(u_file, 'r')
nc_v = Dataset(v_file, 'r')


#Read in file attributes and extract times (Feb 19,2019)

lat = nc_hgt.variables['lat'][:]
lon = nc_hgt.variables['lon'][:] 

narr = {}
time = nc_hgt.variables['time'][:]
timeUnits = nc_hgt.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])

#Assign time index for our time of interest

time_index = np.where(narr['day']==19)[0]

#Now read in meteorological data for this day
level = nc_hgt.variables['level'][:] #In hPa
z = nc_hgt.variables['hgt'][time_index,:,:,:] #In meters
temp = nc_temp.variables['air'][time_index,:,:,:] #In Kelvin
q = nc_q.variables['shum'][time_index,:,:,:]
u = nc_u.variables['uwnd'][time_index,:,:,:]
v = nc_v.variables['vwnd'][time_index,:,:,:]


#%%

#Assign proper units to each variables

temp = temp* units.kelvin
level = level*units.hectopascal







#Define isentropic levels

isentlevs = [296.] * units.kelvin

#Now convert to isentropic coordinates

isentropic = mpcalc.isentropic_interpolation(isentlevs, level, temp, q, u, v, z, tmpk_out = True)


