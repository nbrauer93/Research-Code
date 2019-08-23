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

file = 'GPMCOR_KAR_1407042205_2337_001980_1BS_DAB_05A.h5'


f = h5py.File(file, 'r')
for key in f.keys():
    print(key) #MS HS
    
group = f[key]   


#Print keys inside group
for key in group.keys():
    print(key)
    
#Extract data

lat = group['Latitude'].value     
lon = group['Longitude'].value





    
    
    
    
    
    




 

