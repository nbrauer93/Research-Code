#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 15:42:45 2018

@author: noahbrauer
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
import pickle
import glob
import sys,getopt
from matplotlib.colors import LinearSegmentedColormap

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
    cc = id.variables['CorrelationCoefficient']
    wcc  = id.variables['wCorrelationCoefficient']





	# Create arrays to store binned values
    values    = np.zeros(x['n']*y['n']*z['n'])
    wvalues   = np.zeros(x['n']*y['n']*z['n'])
    kvalues   = np.zeros(x['n']*y['n']*z['n'])
    kwvalues  = np.zeros(x['n']*y['n']*z['n'])
    zvalues   = np.zeros(x['n']*y['n']*z['n'])
    zwvalues  = np.zeros(x['n']*y['n']*z['n'])
    cvalues = np.zeros(x['n']*y['n']*z['n'])
    cwvalues = np.zeros(x['n']*y['n']*z['n'])
    values[:] = float('nan')

	# Add values to arrays
    values[index[:]]  =  (Z_H)[:]
    wvalues[index[:]] = (wZ_H)[:]
    kvalues[index[:]] = (kdp)[:]
    kwvalues[index[:]] = (wkdp)[:]
    zvalues[index[:]] = (zdr)[:]
    zwvalues[index[:]] = (wzdr)[:]
    cvalues[index[:]] = (cc)[:]
    cwvalues[index[:]] = (wcc)[:]



	# Reshape arrays to 3-D GridRad domain
    values  =  values.reshape((z['n'], y['n'] ,x['n']))
    wvalues = wvalues.reshape((z['n'], y['n'] ,x['n']))
    kvalues = kvalues.reshape((z['n'], y['n'] ,x['n']))
    kwvalues = kwvalues.reshape((z['n'], y['n'] ,x['n']))
    zvalues = zvalues.reshape((z['n'], y['n'] ,x['n']))
    zwvalues = zwvalues.reshape((z['n'], y['n'] ,x['n']))
    cvalues = cvalues.reshape((z['n'], y['n'] ,x['n']))
    cwvalues = cwvalues.reshape((z['n'], y['n'] ,x['n']))

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

    cc = {'values'     : cvalues,              \
			 'long_name'  : str(cc.long_name),  \
			 'units'      : str(cc.units),      \
			 'missing'    : float('nan'),        \
			 'wvalues'    : cwvalues,             \
			 'wlong_name' : str(wcc.long_name), \
			 'wunits'     : str(wcc.units),     \
			 'wmissing'   : wcc.missing_value,  \
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
           'cc'                      : cc,   \
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
    
    
file = 'nexrad_3d_v4_0_20170828T211500Z.nc'

nc = Dataset(file ,'r')

lat = nc.variables['Latitude'][:]
lon = nc.variables['Longitude'][:] - 360
alt = nc.variables['Altitude'][:]

lat2,lon2 = np.meshgrid(lat,lon)[:]



files2 = glob.glob('nexrad*.nc')[:]

zdr = np.ones((1441,528,672))*np.nan
kdp = np.ones((1441,528,672))*np.nan
cc = np.ones((1441,528,672))*np.nan
zh = np.ones((1441,528,672))*np.nan

###2 km radar fields

for n,j in enumerate(files2):
    print(j)
    data = read_file(j)

    zdr[n,:,:] = data['zdr']['values'][2,:,:]
    kdp[n,:,:] = data['kdp']['values'][2,:,:]
    zh[n,:,:] = data['Z_H']['values'][2,:,:]
    cc[n,:,:] = data['cc']['values'][2,:,:]    
    
    
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
    diff_reflect_coltbl = LinearSegmentedColormap('DIFF_REFLECT_COLTBL',diff_reflect_cdict)
    return diff_reflect_coltbl
zdrcolor = diff_reflect()




def reflect_ncdc():
    reflect_ncdc_cdict ={
            'red':((0.0000, 0.600, 0.600),
                  (0.0350, 0.450, 0.450),
				(0.0714, 0.000, 0.000),
				(0.1429, 0.000, 0.000),
				(0.2143, 0.000, 0.000),
				(0.2857, 0.000, 0.000),
				(0.3571, 0.000, 0.000),
				(0.4286, 1.000, 1.000),
				(0.5000, 0.906, 0.906),
				(0.5714, 1.000, 1.000),
				(0.6429, 1.000, 1.000),
				(0.7143, 0.839, 0.839),
				(0.7857, 0.753, 0.753),
				(0.8571, 1.000, 1.000),
				(0.9286, 0.600, 0.600),
				(1.000, 0.923, 0.923)),
		'green':	((0.0000, 0.600, 0.600),
                  (0.0350, 0.450, 0.450),
				(0.0714, 0.627, 0.627),
				(0.1429, 0.000, 0.000),
				(0.2143, 1.000, 1.000),
				(0.2857, 0.784, 0.784),
				(0.3571, 0.565, 0.565),
				(0.4286, 1.000, 1.000),
				(0.5000, 0.753, 0.753),
				(0.5714, 0.565, 0.565),
				(0.6429, 0.000, 0.000),
				(0.7143, 0.000, 0.000),
				(0.7857, 0.000, 0.000),
				(0.8571, 0.000, 0.000),
				(0.9286, 0.333, 0.333),
				(1.000, 0.923, 0.923)),
          'blue': ((0.0000, 0.600, 0.600),
                  (0.0350, 0.700, 0.700),
				(0.0714, 0.965, 0.965),
				(0.1429, 0.965, 0.965),
				(0.2143, 0.000, 0.000),
				(0.2857, 0.000, 0.000),
				(0.3571, 0.000, 0.000),
				(0.4286, 0.000, 0.000),
				(0.5000, 0.000, 0.000),
				(0.5714, 0.000, 0.000),
				(0.6429, 0.000, 0.000),
				(0.7143, 0.000, 0.000),
				(0.7857, 0.000, 0.000),
				(0.8571, 1.000, 1.000),
				(0.9286, 0.788, 0.788),
				(1.000, 0.923, 0.923))}
    reflect_ncdc_coltbl = LinearSegmentedColormap('REFLECT_NCDC_COLTBL',reflect_ncdc_cdict)
    return reflect_ncdc_coltbl


zhcolor = reflect_ncdc()
    
    
#%%    
    
    
    
###Calculate 1-hour spatial means    
    
zh_mean = np.ones((120,528,672))*np.nan
zdr_mean = np.ones((120,528,672))*np.nan
kdp_mean = np.ones((120,528,672))*np.nan



for i in range(zh_mean.shape[0]):
    zh_mean[i,:,:] = np.nanmean(zh[12*i:(12*i)+12,:,:],axis = 0)
    zdr_mean[i,:,:] = np.nanmean(zdr[12*i:(12*i)+12,:,:],axis = 0)
    kdp_mean[i,:,:] = np.nanmean(kdp[12*i:(12*i)+12,:,:], axis = 0)
    

    
    
###Plot each hourly mean and store to a GIF:


day = ['8/26','8/27','8/28','8/29']
hour = ['00Z','01Z','02Z','03Z','04Z','05Z','06Z','07Z','08Z','09Z','10Z','11Z','12Z','13Z','14Z','15Z','16Z','17Z','18Z','19Z','20Z','21Z','22Z','23Z']

number = np.arange(0,97,1,dtype = 'int')
alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']


time = []
for i in range(len(day)):
    for j in range(len(hour)):
        junk = str(day[i]) + ' ' + str(hour[j])
        time.append(junk)




        

#%%



for k in range(95):

    cmin = 0.; cmax = 70.; cint = 5.; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=zhcolor,lut=nlevs)
    #plt.figure()
    plt.figure(figsize=(10,6))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
    xlim = np.array([-97,-94]); ylim = np.array([29,30.5])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
    cs = m.contourf(lon2,lat2,zh_mean[k,:,:].T,clevs,cmap=zhcolor,extend='both') #plot lat, lon, and North Pacific SST Anomalies

    #gist_ncar
    #x2star,y2star = m(obs25['lon'],obs25['lat'])
    #m.plot(x2star,y2star,'g*',markersize=2)
    m.drawcounties()
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
        
    x2star,y2star = m(-95.3698,29.7604)
    m.plot(x2star,y2star,'ro',markersize=7, color = 'k')
    label = 'Houston'
    plt.text(x2star+0.05,y2star+0.05,label)    
        
    plt.title('2 km Mean Reflectivity Factor' + ' '  + str(time[k]) + '-' +  str(time[k+1]),name='Calibri',weight='bold',size=16)
    if len(str(number[k]))==1:
        print(number[k])
        plt.savefig('/Users/noahbrauer/Desktop/Research/Radar Data/all/zh_png/Zh0' + str(number[k])+'.png')
    else:
        plt.savefig('/Users/noahbrauer/Desktop/Research/Radar Data/all/zh_png/Zh' +str(number[k])+'.png')
        
    #plt.show(block=False)     
    


import os
import imageio

png_dir = '/Users/noahbrauer/Desktop/Research/Radar Data/all/zh_png/'
images = []
for file_name in os.listdir(png_dir):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('2kmzh.gif', images, duration = 0.75)   


#Now for Zdr

for k in range(1,95):

    cmin = -0.5; cmax = 2.5; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=zhcolor,lut=nlevs)
    #plt.figure()
    plt.figure(figsize=(10,6))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
    xlim = np.array([-97,-94]); ylim = np.array([29,30.5])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
    m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
    m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
    cs = m.contourf(lon2,lat2,zdr_mean[k,:,:].T,clevs,cmap=zdrcolor,extend='both') #plot lat, lon, and North Pacific SST Anomalies

    #gist_ncar
    #x2star,y2star = m(obs25['lon'],obs25['lat'])
    #m.plot(x2star,y2star,'g*',markersize=2)
    m.drawcounties()
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
        
    x2star,y2star = m(-95.3698,29.7604)
    m.plot(x2star,y2star,'ro',markersize=7, color = 'k')
    label = 'Houston'
    plt.text(x2star+0.05,y2star+0.05,label)    
    

        
    plt.title('2 km Mean Differential Reflectivity' + ' '  + str(time[k]) + '-' +  str(time[k+1]),name='Calibri',weight='bold',size=16)
    if len(str(number[k]))==1:
        plt.savefig('/Users/noahbrauer/Desktop/Research/Radar Data/all/zdr_png/Zdr0' +str(number[k])+'.png')
    else:
        plt.savefig('/Users/noahbrauer/Desktop/Research/Radar Data/all/zdr_png/Zdr' +str(number[k])+'.png')
        
    #plt.show(block=False)     
    


import os
import imageio

png_dir = '/Users/noahbrauer/Desktop/Research/Radar Data/all/zdr_png/'
images = []
for file_name in os.listdir(png_dir):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('2kmzdr.gif', images, duration = 0.75)   



    
    
#%%


#Okay, now choose a point (Houston) and plot a time series of 2km reflectvity

houcc = cc[:,228,366]
non_met_index = np.where((houcc>=0.85))[0]


if non_met_index is True:

    houzh = zh[:,228,366]
    houzdr = zdr[:,228,366]
    houkdp = kdp[:,228,366]





plt.figure(figsize=(10,6))
plt.xticks(np.arange(0, 1441, step=160),time[::12], rotation = 45)
#plt.xlabel('Time', size = 14)
plt.ylabel('dBZ', size = 14)
plt.title('Houston 2 km $Z_{H}$', size = 16)
plt.plot(houzh, color = 'g')
plt.show()



plt.figure(figsize=(10,6))
#plt.xlabel('Time')
plt.ylabel('dB', size = 14)
plt.title('Houston 2 km $Z_{DR}$',size = 16)
plt.xticks(np.arange(0, 1441, step=160),time[::12], rotation = 45)
plt.plot(houzdr, color = 'r')
plt.show()
    
    
plt.figure(figsize=(10,6))
#plt.xlabel('Time', size = 14)
plt.ylabel('deg/km', size = 14)
plt.title('Houston 2 km $K_{DP}$',size = 16)
plt.xticks(np.arange(0, 1441, step=160),time[::12], rotation = 45)
plt.ylim(-1,5)
plt.plot(houkdp, color = 'b')
plt.show()        



   
    
    
    
    
    