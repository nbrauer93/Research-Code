# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 00:25:44 2018

@author: noahb
"""

import matplotlib.pyplot as plt 
from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from mpl_toolkits.basemap import Basemap,maskoceans,interp,shiftgrid
from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator

ncFile = 'sst2018.nc'

narr = {}

with Dataset(ncFile,'r') as nc:
    lon = nc.variables['lon'][:]-180
    lat = nc.variables['lat'][:]
    #lam = nc.variables['Lambert_Conformal'][:]
    time = nc.variables['time'][:]
    anom = nc.variables['anom'][:].T #shape: 349/277/2920 = IJT convert from mm to in
    #narr['lat'],narr['lon']=np.meshgrid(lat,lon)
    timeUnits = nc.variables['time'].units
    tmpDates = num2date(time,timeUnits,calendar='gregorian')
    narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
    narr['day'] = np.asarray([d.day for d in narr['date']])
    narr['month'] = np.asarray([d.month for d in narr['date']]) 
    narr['year'] = np.asarray([d.year for d in narr['date']])
    
    anom2 = np.where(anom.mask,np.nan,anom.data)

aug_index = np.where((narr['month']==5))[0] 
anom_may = anom2[:,:,aug_index] #349, 277, 248

anom_total = np.nanmean(anom_may,axis = 2)


#%%

cmin = 0; cmax = 8.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
#plt.figure()
#plt.figure(figsize=(24,16))
#cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
#plt.figure(figsize=(18,18))
xlim = np.array([-140.,-60.]); ylim = np.array([10.,45.])
#parallels = np.arange(23.,35.,1.)
# labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
#xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,anom_total.T,clevs,cmap='GnBu',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
m.drawcounties()
#x2star,y2star = m(obs25['lon'],obs25['lat'])
#m.plot(x2star,y2star,'g*',markersize=2)

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('in',weight='bold',name='Calibri',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('Sea surface temperature anomalies',name='Calibri',weight='bold',size=16)
'''
x2star,y2star = m(-97.5164,35.4676)
m.plot(x2star,y2star,'ro',markersize=7)
label = 'Oklahoma City'

#plt.text(x2star+0.1,y2star+0.1,label)
label = 'Image by Noah Brauer'
plt.text(-102.5,32.5,label)
label2 = 'Data Provided by NOAA/ESRL'
plt.text(-102.5,32.25,label2)
'''
    

plt.show(block=False) 