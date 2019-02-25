
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np

from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import pickle
import glob
import sys,getopt
from matplotlib.colors import LinearSegmentedColormap


import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap



file = '20170827-040115.netcdf'

nc = Dataset(file, 'r')


lat = np.linspace(55.0005,20.0005,7000)     
lon = np.linspace(-130.0005,-60.0005,14000)
lat2,lon2 = np.meshgrid(lat,lon)


rotation = nc.variables['RotationTrack60min'][:] #Micheal's file 







plt.figure(figsize=(18,12))
xlim = np.array([-97,-94]); ylim = np.array([27,32])
#xlim = np.array([-130.0005,-60.0005]); ylim = np.array([20.0005,55.0005])
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
   
cs = m.contour(lon2.T,lat2.T,rotation, cmap = 'YlOrRd') 
m.drawcounties()

plt.show()
