# -*- coding: utf-8 -*-
"""
Created on Sun May  6 19:44:39 2018

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
from colormap import Colormap
import glob

#from gridrad.py import read_file

def read_file(infile):
	
	# Import python libraries
    import sys
    import os
    import numpy as np
    import netCDF4
	
	# Check to see if file exists
    if not os.path.isfile(infile):
        print('File "' + infile + '" does not exist.  Returning -2.')
        return -2
	
	# Check to see if file has size of zero
    if os.stat(infile).st_size == 0:
        print('File "' + infile + '" contains no valid data.  Returning -1.')
        return -1
	
    from netCDF4 import Dataset
    from netCDF4 import Variable
	# Open GridRad netCDF file
    id = Dataset(infile, "r", format="NETCDF4")
	
	# Read global attributes
    Analysis_time           = str(id.getncattr('Analysis_time'          ))
    Analysis_time_window    = str(id.getncattr('Analysis_time_window'   ))
    File_creation_date      = str(id.getncattr('File_creation_date'     ))
    Grid_scheme             = str(id.getncattr('Grid_scheme'            ))
    Algorithm_version       = str(id.getncattr('Algorithm_version'      ))
    Algorithm_description   = str(id.getncattr('Algorithm_description'  ))
    Data_source             = str(id.getncattr('Data_source'            ))
    Data_source_URL         = str(id.getncattr('Data_source_URL'        ))
    NOAA_wct_export_Version = str(id.getncattr('NOAA_wct-export_Version'))
    Authors                 = str(id.getncattr('Authors'                ))
    Project_sponsor         = str(id.getncattr('Project_sponsor'        ))
    Project_name            = str(id.getncattr('Project_name'           ))
	
	# Read list of merged files
    file_list    = (id.variables['files_merged'])[:]
    files_merged = ['']*(id.dimensions['File'].size)
    for i in range(0,id.dimensions['File'].size):
        for j in range(0,id.dimensions['FileRef'].size):
            files_merged[i] += str(file_list[i,j])
	
	# Read longitude dimension
    x = id.variables['Longitude']
    x = {'values'    : x[:],             \
		  'long_name' : str(x.long_name), \
		  'units'     : str(x.units),     \
		  'delta'     : str(x.delta),     \
		  'n'         : len(x[:])}
	
	# Read latitude dimension
    y = id.variables['Latitude']
    y = {'values'    : y[:],             \
		  'long_name' : str(y.long_name), \
		  'units'     : str(y.units),     \
		  'delta'     : str(y.delta),     \
		  'n'         : len(y[:])}
	
	# Read altitude dimension
    z = id.variables['Altitude']
    z = {'values'    : z[:],             \
		  'long_name' : str(z.long_name), \
		  'units'     : str(z.units),     \
		  'delta'     : str(z.delta),     \
		  'n'         : len(z[:])}
	
	# Read observation and echo counts
    nobs  = (id.variables['Nradobs' ])[:]
    necho = (id.variables['Nradecho'])[:]
    index = (id.variables['index'   ])[:]
	
	# Read reflectivity variables	
    Z_H  = id.variables['Reflectivity' ]
    wZ_H = id.variables['wReflectivity']
    zdr = id.variables['DifferentialReflectivity']
    kdp = id.variables['DifferentialPhase']
    wzdr = id.variables['wDifferentialReflectivity']
    wkdp = id.variables['wDifferentialPhase']
    
    
    
     

	# Create arrays to store binned values	
    values    = np.zeros(x['n']*y['n']*z['n'])
    wvalues   = np.zeros(x['n']*y['n']*z['n'])
    kvalues   = np.zeros(x['n']*y['n']*z['n'])
    kwvalues  = np.zeros(x['n']*y['n']*z['n'])
    zvalues   = np.zeros(x['n']*y['n']*z['n'])
    zwvalues  = np.zeros(x['n']*y['n']*z['n'])
    values[:] = float('nan')

	# Add values to arrays
    values[index[:]]  =  (Z_H)[:]
    wvalues[index[:]] = (wZ_H)[:]
    kvalues[index[:]] = (kdp)[:]
    kwvalues[index[:]] = (wkdp)[:]
    zvalues[index[:]] = (zdr)[:]
    zwvalues[index[:]] = (wzdr)[:]
    
    
	
	# Reshape arrays to 3-D GridRad domain
    values  =  values.reshape((z['n'], y['n'] ,x['n']))
    wvalues = wvalues.reshape((z['n'], y['n'] ,x['n']))
    kvalues = kvalues.reshape((z['n'], y['n'] ,x['n']))
    kwvalues = kwvalues.reshape((z['n'], y['n'] ,x['n']))
    zvalues = zvalues.reshape((z['n'], y['n'] ,x['n']))
    zwvalues = zwvalues.reshape((z['n'], y['n'] ,x['n']))
    

    Z_H = {'values'     : values,              \
			 'long_name'  : str(Z_H.long_name),  \
			 'units'      : str(Z_H.units),      \
			 'missing'    : float('nan'),        \
			 'wvalues'    : wvalues,             \
			 'wlong_name' : str(wZ_H.long_name), \
			 'wunits'     : str(wZ_H.units),     \
			 'wmissing'   : wZ_H.missing_value,  \
			 'n'          : values.size}
    
    zdr = {'values'     : zvalues,              \
			 'long_name'  : str(zdr.long_name),  \
			 'units'      : str(zdr.units),      \
			 'missing'    : float('nan'),        \
			 'wvalues'    : zwvalues,             \
			 'wlong_name' : str(wzdr.long_name), \
			 'wunits'     : str(wzdr.units),     \
			 'wmissing'   : wzdr.missing_value,  \
			 'n'          : values.size}
    
    kdp = {'values'     : kvalues,              \
			 'long_name'  : str(kdp.long_name),  \
			 'units'      : str(kdp.units),      \
			 'missing'    : float('nan'),        \
			 'wvalues'    : kwvalues,             \
			 'wlong_name' : str(wkdp.long_name), \
			 'wunits'     : str(wkdp.units),     \
			 'wmissing'   : wkdp.missing_value,  \
			 'n'          : values.size}
	
	# Close netCDF4 file
    id.close()
	
	# Return data dictionary	
    return {'name'                    : 'GridRad analysis for ' + Analysis_time, \
			  'x'                       : x, \
			  'y'                       : y, \
			  'z'                       : z, \
			  'Z_H'                     : Z_H, \
           'zdr'                     : zdr, \
           'kdp'                     : kdp, \
			  'nobs'                    : nobs, \
			  'necho'                   : necho, \
			  'file'                    : infile, \
			  'files_merged'            : files_merged, \
			  'Analysis_time'           : Analysis_time, \
			  'Analysis_time_window'    : Analysis_time_window, \
			  'File_creation_date'      : File_creation_date, \
			  'Grid_scheme'             : Grid_scheme, \
			  'Algorithm_version'       : Algorithm_version, \
			  'Algorithm_description'   : Algorithm_description, \
			  'Data_source'             : Data_source, \
			  'Data_source_URL'         : Data_source_URL, \
			  'NOAA_wct_export_Version' : NOAA_wct_export_Version, \
			  'Authors'                 : Authors, \
			  'Project_sponsor'         : Project_sponsor, \
			  'Project_name'            : Project_name}
    
    

#%%

files2 = glob.glob('nexrad*.nc')[:]

zdr = np.ones((145,528//3,672//3))*np.nan
kdp = np.ones((145,528//3,672//3))*np.nan
zh = np.ones((145,528//3,672//3))*np.nan


for n,j in enumerate(files2):
    data = read_file(j)

    zdr[n,:,:] = data['zdr']['values'][4,::3,::3]
    kdp[n,:,:] = data['kdp']['values'][4,::3,::3]
    zh[n,:,:] = data['Z_H']['values'][4,::3,::3]
    

    
#%%
    
ncFile = 'apcp.2017.nc'

narr = {}

with Dataset(ncFile,'r') as nc:
    lon = nc.variables['lon'][:]
    #lon2 = lon+360
    lat = nc.variables['lat'][:]
    lam = nc.variables['Lambert_Conformal'][:]
    time = nc.variables['time'][:]
    precip = nc.variables['apcp'][:].T #shape: 349/277/2920 = IJT
    #narr['lat'],narr['lon']=np.meshgrid(lat,lon)
    timeUnits = nc.variables['time'].units
    tmpDates = num2date(time,timeUnits,calendar='gregorian')
    narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
    narr['day'] = np.asarray([d.day for d in narr['date']])
    narr['month'] = np.asarray([d.month for d in narr['date']])
    narr['year'] = np.asarray([d.year for d in narr['date']])
    
    precip2 = np.where(precip.mask,np.nan,precip.data)



aug_index = np.where((narr['month']==8) & (narr['day']>25) & (narr['day']<28))[0]
precip_aug = precip2[:,:,aug_index] #349, 277, 248

I,J,T = precip_aug.shape
precip_squish = precip_aug.reshape(I*J,T,order='F') #(96673, 16)   
precip = precip_squish.T 

#Standardize precip
precip_ano = precip-np.nanmean(precip,0)
precip = precip_ano/np.nanstd(precip,0)

nanremove = np.where(~np.isnan(precip[0,:]))[0]

precipraw = precip[:,nanremove] ###raw precip with no nans (not anomalies)

precip3 = precip[:,nanremove] #Goes into EOF machine
precip_ano_no_nan = precip_ano[:,nanremove] #Goes into regression

#%%
preciptotal = np.zeros((2,96673))*np.nan

for i in range(0,2):
    preciptotal[i,:] = np.sum(precip[i*8:8*i+8,:],axis = 0)

T2 = 2
I2 = 277
J2 = 349
preciptotal2 = preciptotal.reshape(T2,I2,J2)



  
    
#%%    
    
    
    
cmin = 0; cmax = 10.; cint = 0.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
#plt.figure()
#plt.figure(figsize=(24,16))
#cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-101,-92]); ylim = np.array([26,33])
#parallels = np.arange(23.,35.,1.)
# labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates() #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
#xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon,lat,preciptotal2[0],clevs,cmap='GnBu',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
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
plt.title('Total Precipitation (in) 2017082600Z-2017082800Z',name='Calibri',weight='bold',size=16)
x2star,y2star = m(-97.3964,27.8006)
m.plot(x2star,y2star,'ro',markersize=7, label = 'Corpus Christi')
    
x3star,y3star = m(-95.3698,29.7604)
m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
plt.legend(loc = 4)
plt.show(block=False) 





#%%
###Take three hour averages of zdr and kdp (16 3-hour averages)

#zdravg = np.zeros((16,354816))
#kdpavg = np.zeros((16,354816))

zdravg = np.ones_like(zdr)[:16,:,:]
kdpavg = np.ones_like(kdp)[:16,:,:]

for i in range(0,16):
    zdravg[i,:,:] = np.nanmean(zdr[i*9:9*i+9,:],0)
    kdpavg[i,:,:] = np.nanmean(kdp[i*9:9*i+9,:],0)

###Standardize the data

zdr3 = (zdravg-np.nanmean(zdravg,0))/np.nanstd(zdravg,0)
kdp3 = (kdpavg-np.nanmean(kdpavg,0))/np.nanstd(kdpavg,0) #(16,354816)

T,J,I = zdravg.shape


#Import the lat-lon for radar

radarlon = data['x']['values'][::3]-360
radarlat = data['y']['values'][::3]

radarlat2,radarlon2 = np.meshgrid(radarlat,radarlon)

zh1 = zh[0,:,:] #8/26/00Z
zh2 = zh[72,:,:] #8/27/00Z
zh3 = zh[-1,:,:] #8/28/00Z

zdralone1 = zdr[0,:,:]
zdralone2 = zdr[72,:,:]
zdralone3 = zdr[-1,:,:]





#%%

cmin = 0; cmax = 75.; cint = 2; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
#plt.figure()
#plt.figure(figsize=(24,16))
#cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-101,-92]); ylim = np.array([26,33])
#parallels = np.arange(23.,35.,1.)
# labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates() #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
#xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(radarlon2,radarlat2,zh1.T,clevs,cmap='gist_ncar',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
#x2star,y2star = m(obs25['lon'],obs25['lat'])
#m.plot(x2star,y2star,'g*',markersize=2)

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dBZ',weight='bold',name='Calibri',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('3 km Reflectivity (dBZ) 2017082600Z',name='Calibri',weight='bold',size=16)
x2star,y2star = m(-97.3964,27.8006)
m.plot(x2star,y2star,'r*',markersize=10, label = 'Corpus Christi')
    
x3star,y3star = m(-95.3698,29.7604)
m.plot(x3star,y3star,'k*',markersize=10, label = 'Houston')
plt.legend(loc = 4)
plt.show(block=False) 

#%%


def diff_reflect():
    diff_reflect_cdict ={
            'red':((0.000, 0.000, 0.000),
                   (0.333, 1.000, 1.000),
				 (0.417, 0.000, 0.000),
                   (0.500, 0.000, 0.000),
                   (0.583, 1.000, 1.000),
                   (0.750, 1.000, 1.000),
                   (0.833, 1.000, 1.000),
                   (1.000, 1.000, 1.000)),
		'green':	((0.000, 0.000, 0.000),
                   (0.333, 1.000, 1.000),
				 (0.417, 0.000, 0.000),
                   (0.500, 1.000, 1.000),
                   (0.583, 1.000, 1.000),
                   (0.750, 0.000, 0.000),
                   (0.833, 0.000, 0.000),
                   (1.000, 1.000, 1.000)),
          'blue': ((0.000, 0.000, 0.000),
                   (0.333, 1.000, 1.000),
				 (0.417, 1.000, 1.000),
                   (0.500, 1.000, 1.000),
                   (0.583, 0.000, 0.000),
                   (0.750, 0.000, 0.000),
                   (0.833, 1.000, 1.000),
                   (1.000, 1.000, 1.000))}
    diff_reflect_coltbl = Colormap('DIFF_REFLECT_COLTBL',diff_reflect_cdict)  #You will need to import this
    return diff_reflect_coltbl

zdrcmap = diff_reflect


#%%

cmin = -6; cmax = 6.; cint = 0.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
#plt.figure()
#plt.figure(figsize=(24,16))
#cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-101,-92]); ylim = np.array([26,33])
#parallels = np.arange(23.,35.,1.)
# labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates() #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
#xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(radarlon2,radarlat2,zdralone1.T,clevs,cmap='gist_ncar',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
#x2star,y2star = m(obs25['lon'],obs25['lat'])
#m.plot(x2star,y2star,'g*',markersize=2)

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dB',weight='bold',name='Calibri',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('3 km $Z_{DR}$ (dB) 2017082600Z',name='Calibri',weight='bold',size=16)
x2star,y2star = m(-97.3964,27.8006)
m.plot(x2star,y2star,'ro',markersize=7, label = 'Corpus Christi')
    
x3star,y3star = m(-95.3698,29.7604)
m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
plt.legend(loc = 4)
plt.show(block=False) 







#%%


T,I,J = zdravg.shape
zdravg2 = zdravg.reshape((T,I*J), order = 'F')
kdpavg2 = kdpavg.reshape((T,I*J), order = 'F')


newzdr = np.zeros((16,39424))  
newkdp = np.zeros((16,39424))


for n in range(0,16):
    for i in range(0,39424):
        if precip_ano[n,i]>=0.05:
            newzdr[n,i] = zdravg2[n,i]
            newkdp[n,i] = kdpavg2[n,i]
            
        else:
            newzdr[n,i] = np.nan
            newkdp[n,i] = np.nan
            
zdrmean = np.nanmean(newzdr, axis=0)
kdpmean = np.nanmean(newkdp, axis=0)        
        

I2 = 176
J2 = 224
zdrmean1 = zdrmean.reshape((I2,J2), order = 'F')
kdpmean1 = kdpmean.reshape((I2,J2), order = 'F')



cmin = -5; cmax = 5.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
#plt.figure()
#plt.figure(figsize=(24,16))
#cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-101,-92]); ylim = np.array([26,33])

#parallels = np.arange(23.,35.,1.)
# labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates() #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])


#xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(radarlon2,radarlat2,zdrmean1.T,clevs,cmap='RdPu',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
#x2star,y2star = m(obs25['lon'],obs25['lat'])
#m.plot(x2star,y2star,'g*',markersize=2)

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('dB',weight='bold',name='Calibri',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('3 km $Z_{DR}$ (dB) (Precip anomalies > 0.05")',name='Calibri',weight='bold',size=16)
x2star,y2star = m(-97.3964,27.8006)
m.plot(x2star,y2star,'bo',markersize=7, label = 'Corpus Christi')
    
x3star,y3star = m(-95.3698,29.7604)
m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
plt.legend(loc = 4)
plt.show(block=False) 



cmin = -2.; cmax = 2.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
#plt.figure()
#plt.figure(figsize=(24,16))
#cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-101,-92]); ylim = np.array([26,33])

#parallels = np.arange(23.,35.,1.)
# labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries(), m.drawstates() #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])


#xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(radarlon2,radarlat2,kdpmean1.T,clevs,cmap='Greens',extend='both') #plot lat, lon, and North Pacific SST Anomalies 
#x2star,y2star = m(obs25['lon'],obs25['lat'])
#m.plot(x2star,y2star,'g*',markersize=2)

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('deg/km',weight='bold',name='Calibri',size=14)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
cbar.set_ticks(clevs[::4])
cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(14)
plt.title('3 km $K_{DP}$ (deg/km) (Precip anomalies > 0.05")',name='Calibri',weight='bold',size=16)
x2star,y2star = m(-97.3964,27.8006)
m.plot(x2star,y2star,'ro',markersize=7, label = 'Corpus Christi')
    
x3star,y3star = m(-95.3698,29.7604)
m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
plt.legend(loc = 4)
plt.show(block=False) 
       
        
        
###Test for significance

zdrmean2 = np.nanmean(newzdr)
kdpmean2 = np.nanmean(newkdp)


print('The mean Kdp value when precipitation anomalies exceed 0.05" in 3-hours is',kdpmean2,'deg/km')        
print('The mean Zdr value when precipitation anomalies exceed 0.05" in 3-hours is',zdrmean2,'dB')        
        
zdrvariance = (np.nanstd(newzdr))**2
kdpvariance = (np.nanstd(newkdp))**2

print('The variance of Kdp when precipitation anomalies exceed 0.05" in 3-hours is',kdpvariance,'deg/km')        
print('The variance of Zdr when precipitation anomalies exceed 0.05" in 3-hours is',zdrvariance,'dB')        
        


###two tail t-test


##Null hypothesis: 


T,I,J = zdr3.shape
zdrcomp = zdr3.reshape(T,I*J) #zdr for all precip

zdrsig = scipy.stats.ttest_ind(zdrcomp,newzdr, nan_policy='omit', axis = 0)[1] 
print(zdrsig)
        
        
        
        
        
        
        









