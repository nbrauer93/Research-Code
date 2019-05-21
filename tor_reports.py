#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:41:31 2019

@author: noahbrauer
"""

import numpy as np
import matplotlib.pyplot as plt

from numpy import genfromtxt


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap


file = '2017_torn.csv'

my_data = genfromtxt(file, delimiter=',')

reports = {}

reports['month'] = my_data[:,2]
reports['day'] = my_data[:,3]
reports['lat'] = my_data[:,17]
reports['lon'] = my_data[:,18]


august_index = np.where((reports['month']==8)&(reports['day']>25))[0]

lat = reports['lat'][august_index]
lon = reports['lon'][august_index]

coords = np.array([lat,lon], dtype = 'float').T

lat2,lon2 = np.meshgrid(lat,lon)[:]


###Now we plot

plt.figure(figsize=(10,10))
xlim = np.array([-98,-87]); ylim = np.array([25,37])


m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
m.drawcounties()
cs = m(lon2,lat2) 
plt.title('Filtered Tornado Reports 26-31 August 2017', size = 22, weight = 'bold', name = 'Calibri')

tor_x, tor_y = m(coords[:,1], coords[:,0])
m.scatter(tor_x,tor_y, color = 'r')







plt.show()

