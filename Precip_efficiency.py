#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:12:16 2018

@author: noahbrauer
"""

from ecmwfapi import ECMWFDataServer
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.interpolate import RectBivariateSpline
from netCDF4 import Dataset, num2date, MFDataset
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
import glob
from scipy.interpolate import griddata
from numpy.random import uniform, seed



import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

file = 'flux.nc'
file2 = 'precip.nc'
file3 = 'specific_hum.nc'

era = {}

nc = Dataset(file,'r')
nc2 = Dataset(file2,'r')
nc3 = Dataset(file3,'r')

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


moisture_flux_div = nc.variables['p84.162'][:] #units of kgm^-2s^-1
total_precip = nc2.variables['tp'][:] #units of m
q = nc3.variables['q'][:]


aug_index = np.where(era['day']>=24)[0]

moisture_flux_div_aug = moisture_flux_div[aug_index,:,:]
total_precip_aug = total_precip[aug_index,:,:]
q_aug = q[aug_index,:,:]
const = -1/9.8 #1/g


#Now read in the PRISM Data


prism_path = glob.glob('*_bil.bil')

 
def read_prism_hdr(hdr_path):
    """Read an ESRI BIL HDR file"""
    with open(hdr_path, 'r') as input_f:
        header_list = input_f.readlines()
    return dict(item.strip().split() for item in header_list)
 
def read_prism_bil(bil_path):
    """Read an array from ESRI BIL raster file"""
    hdr_dict = read_prism_hdr(bil_path.replace('.bil', '.hdr'))
    # For now, only use NROWS, NCOLS, and NODATA
    # Eventually use NBANDS, BYTEORDER, LAYOUT, PIXELTYPE, NBITS
 
    prism_array = np.fromfile(bil_path, dtype=np.float32)
    prism_array = prism_array.reshape(
        int(hdr_dict['NROWS']), int(hdr_dict['NCOLS']))
    prism_array[prism_array == float(hdr_dict['NODATA'])] = np.nan
    return prism_array
 
#prism_path = 'PRISM_ppt_stable_4kmD2_20170828_bil.bil'
    
prism_array = np.ones((621,1405,5))*np.nan

for i in range(len(prism_path)):
    prism_array[:,:,i] = read_prism_bil(prism_path[i])

prism_array = prism_array.T


# making the data grid points for PRISM

xlim = np.array([-125.0208333, -66.4791667])
ylim = np.array([24.0625000, 49.9375000])
prism_cols = 1405
prism_rows = 621

x = np.linspace(-125.0208333, -66.4791667, 1405)
y = np.linspace(24.0625000, 49.9375000, 621)[::-1]
 
iX = np.where( (x>=xlim[0]) & (x<=xlim[1]) )[0]
iY = np.where( (y>=ylim[0]) & (y<=ylim[1]) )[0]
prism_nodata = -9999
lat3,lon3 = np.meshgrid(y[iY],x[iX])


#total_precip  = np.nansum(prism_array, axis = 2)

#%%
##Take derivatives of q: 
                                                                                                                                                                                                                 

#Now integrate for each time (from 1000 to 300 mb; sum up all values)


delta_q = np.ones((7,20,241,480))*np.nan
for i in range(q_aug.shape[0]-1):
    delta_q[i,:,:,:] = q_aug[i+1,:,:,:] - q_aug[i,:,:,:]

#Sum the vertical dimension (integrate from 1000 mb to 300 mb)
        

delta_q_int = np.sum(delta_q, axis = 1)

#Now calculate precipitation efficiency:
 
bottom_sum = const*(delta_q_int + moisture_flux_div_aug[1:,:,:])


#%%
###Now transform PRISM data to the coarser, ERA-interim data grid


num_prism_times = prism_array.shape[0]
prism_array_interp = np.full((num_prism_times, len(lon), len(lat)), np.nan)

for i in range(num_prism_times):
    this_precip_matrix = np.fliplr(prism_array[i, ...])
    # this_precip_matrix[np.isnan(this_precip_matrix)] = 0.
    
    interp_object = RectBivariateSpline(x, y[::-1], this_precip_matrix, kx=1, ky=1, s=0)
    prism_array_interp[i, ...] = interp_object(lon, lat[::-1], grid=True)
    prism_array_interp[i, ...] = np.fliplr(prism_array_interp[i, ...])


T,I,J = prism_array_interp.shape

prism_array_interp_reshape = prism_array_interp.reshape(T,J,I, order = 'F')/1000  #convert from mm to m



#precip_eff = (prism_array_interp_reshape[:,:,:]/bottom_sum[1:-1,:,:])
'''

precip_eff = np.ones((6,241,480))*np.nan

for k in range(precip_eff.shape[0]):
    precip_eff[k,:,:] = (total_precip_aug[1+k:k-1,:,:]/bottom_sum[1+k:,:,:])  ###25-30
'''    



   
    
    
    


#%%

'''
#Plot

plt.figure(figsize=(10,6))
cmin = 0; cmax = 100.; cint = 2.; clevs = np.arange(cmin,cmax,cint)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
xlim = np.array([-99.,-93.]); ylim = np.array([27.,33.])
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates(),  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
cs = m.contourf(lon2,lat2,precip_eff[4,:,:].T,cmap='BrBG',extend='both') 
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('%',weight='bold',name='Calibri',size=14)
cticks = []
plt.title('Precipitation Efficiency',name='Calibri',weight='bold',size=16)
plt.show(block=False)
'''


#%%


#Okay, now let's bring in the NARR files. 

precip = 'apcp.2017.nc'
pwat = 'pr_wtr.2017.nc'
flux = 'mconv.hl1.2017.nc'


pwat_load = Dataset(pwat, 'r')
flux_load = Dataset(flux,'r')


narr = {}

with Dataset(precip,'r') as ncnarr:
    lon = ncnarr.variables['lon'][:]
    #lon2 = lon+360
    lat = ncnarr.variables['lat'][:]
    lam = ncnarr.variables['Lambert_Conformal'][:]
    time = ncnarr.variables['time'][:]
    precip = ncnarr.variables['apcp'][:] #shape: 349/277/2920 = IJT convert from mm to in
    #narr['lat'],narr['lon']=np.meshgrid(lat,lon)
    timeUnits = ncnarr.variables['time'].units
    tmpDates = num2date(time,timeUnits,calendar='gregorian')
    narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
    narr['day'] = np.asarray([d.day for d in narr['date']])
    narr['month'] = np.asarray([d.month for d in narr['date']])
    narr['year'] = np.asarray([d.year for d in narr['date']])
    
    precip2 = np.where(precip.mask,np.nan,precip.data)

pwat_data = pwat_load.variables['pr_wtr'][:]
mflux = flux_load.variables['mconv'][:]


august_narr_index = np.where((narr['month']==8)&(narr['day']>24))[0]

pwat_aug = pwat_data[august_narr_index,:,:]
precip2_aug = precip2[august_narr_index,:,:] 
mflux_aug = mflux[august_narr_index,:,:]
 


efficiency = (precip2_aug/pwat_aug)*100

#%%

cmin = 0.; cmax = 400.; cint = 25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,6))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-99,-92]); ylim = np.array([25,32])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,efficiency[5,:,:],clevs,cmap='Greens',extend='both') #plot lat, lon, and North Pacific SST Anomalies

m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('%',weight='bold',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)

plt.title('Precipitation Efficiency 8/30',name='Calibri',size=16, weight = 'bold')
plt.show(block=False)

#%%
###Correlate precip efficiency and moisture flux convergence (mflux_aug, efficiency)


corr = np.ones((7,277,349))*np.nan

for i in range(0,277):
    for j in range(0,349):
        corr[:,i,j] = stats.pearsonr(efficiency[:,i,j],mflux_aug[:,i,j])[0]
#%%

#Okay, now read in wind data for streamlines
        
wind_u = 'uwnd.10m.2017.nc'  
wind_v = 'vwnd.10m.2017.nc'

ncu = Dataset(wind_u,'r')
ncv = Dataset(wind_v,'r')

u = ncu.variables['uwnd'][:]
v = ncv.variables['vwnd'][:]      
        
u_aug = u[august_narr_index,:,:]
v_aug = v[august_narr_index,:,:]
'''
lat1 = lat[:,0]
lon1 = lon[0,:]

lons,lats = np.meshgrid(lat1,lon1)
'''



        
#%%
#Plot correlation


cmin = -1.; cmax = 1.1; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(10,6))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-99,-92]); ylim = np.array([25,32])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,corr[5,:,:],clevs,cmap='coolwarm',extend='both') #plot lat, lon, and North Pacific SST Anomalies
m.drawcounties()
m.quiver(lon,lat,u_aug[5,:,:], v_aug[5,:,:],scale = 400)
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('Correlation',weight='bold',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)


    


plt.title('Precipitation Efficiency vs. Horizontal Moisture Convergence 8/29',name='Calibri',size=14, weight = 'bold')
plt.show(block=False)


