#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:25:07 2018

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
    
file = 'nexrad_3d_v4_0_20170827T000000Z.nc'

nc = Dataset(file ,'r')

lat = nc.variables['Latitude'][:]
lon = nc.variables['Longitude'][:] - 360
alt = nc.variables['Altitude'][:]


####Houston point 227,365 (lat,lon)


files2 = glob.glob('nexrad*.nc')[:]

zdr = np.ones((5,28))*np.nan
kdp = np.ones((5,28))*np.nan
cc = np.ones((5,28))*np.nan
zh = np.ones((5,28))*np.nan

###3 km radar fields

for n,j in enumerate(files2):
    print(j)
    data = read_file(j)

    zdr[n,:,:] = data['zdr']['values'][:,227,365]
    kdp[n,:,:] = data['kdp']['values'][:,227,365]
    zh[n,:,:] = data['Z_H']['values'][:,227,365]
    cc[n,:,:] = data['cc']['values'][:,227,365]     



#%%


















####Rhohv

fig2 = plt.figure(figsize = (18,12))
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

ax1.plot(cc2,alt, color = 'g')
ax2.plot(cc3, alt, color = 'g')
ax3.plot(cc4, alt, color = 'g')
ax4.plot(cc5, alt, color = 'g')
ax5.plot(cc6, alt, color = 'g')

ax1.set_ylim(1,10)
ax2.set_ylim(1,10)
ax3.set_ylim(1,10)
ax4.set_ylim(1,10)
ax5.set_ylim(1,10)

ax1.set_xlim(0.8,1)
ax2.set_xlim(0.8,1)
ax3.set_xlim(0.8,1)
ax4.set_xlim(0.8,1)
ax5.set_xlim(0.8,1)


ax1.set_xlabel(r'$\rho_{HV}$', size = 12)
ax2.set_xlabel(r'$\rho_{HV}$',size = 12)
ax3.set_xlabel(r'$\rho_{HV}$',size = 12)
ax4.set_xlabel(r'$\rho_{HV}$',size = 12)
ax5.set_xlabel(r'$\rho_{HV}$',size = 12)

ax1.set_ylabel('Altitude (km)', size = 12)
ax2.set_ylabel('Altitude (km)',size = 12)
ax3.set_ylabel('Altitude (km)',size = 12)
ax4.set_ylabel('Altitude (km)',size = 12)
ax5.set_ylabel('Altitude (km)',size = 12)


ax1.set_title(' 826 00Z', size = 16)
ax2.set_title('8/26 12Z', size = 16)
ax3.set_title(' 8/27 00Z', size = 16)
ax4.set_title('8/27 12Z', size = 16)
ax5.set_title('8/28 00Z', size = 16)
subplots_adjust(wspace = 0.4, hspace = 0.3)

fig2 = gcf()
st = fig2.suptitle(r"Houston $\rho_{HV}$", size=24)
plt.subplots()


###ZDR

fig2 = plt.figure(figsize = (18,12))
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

ax1.plot(zdr2,alt, color = 'r')
ax2.plot(zdr3, alt, color = 'r')
ax3.plot(zdr4, alt, color = 'r')
ax4.plot(zdr5, alt, color = 'r')
ax5.plot(zdr6, alt, color = 'r')

ax1.set_ylim(1,10)
ax2.set_ylim(1,10)
ax3.set_ylim(1,10)
ax4.set_ylim(1,10)
ax5.set_ylim(1,10)

ax1.set_xlabel('dB', size = 12)
ax2.set_xlabel('dB',size = 12)
ax3.set_xlabel('dB',size = 12)
ax4.set_xlabel('dB',size = 12)
ax5.set_xlabel('dB',size = 12)

ax1.set_ylabel('Altitude (km)', size = 12)
ax2.set_ylabel('Altitude (km)',size = 12)
ax3.set_ylabel('Altitude (km)',size = 12)
ax4.set_ylabel('Altitude (km)',size = 12)
ax5.set_ylabel('Altitude (km)',size = 12)


ax1.set_title(' 826 00Z', size = 16)
ax2.set_title('8/26 12Z', size = 16)
ax3.set_title(' 8/27 00Z', size = 16)
ax4.set_title('8/27 12Z', size = 16)
ax5.set_title('8/28 00Z', size = 16)
subplots_adjust(wspace = 0.4, hspace = 0.3)

fig2 = gcf()
st = fig2.suptitle("Houston $Z_{DR}$", size=24)
plt.subplots()









###KDP
fig2 = plt.figure(figsize = (18,12))
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

ax1.plot(kdp2,alt, color = 'k')
ax2.plot(kdp3, alt, color = 'k')
ax3.plot(kdp4, alt, color = 'k')
ax4.plot(kdp5, alt, color = 'k')
ax5.plot(kdp6, alt, color = 'k')

ax1.set_ylim(1,10)
ax2.set_ylim(1,10)
ax3.set_ylim(1,10)
ax4.set_ylim(1,10)
ax5.set_ylim(1,10)

ax1.set_xlim(-1,5)
ax2.set_xlim(-1,5)
ax3.set_xlim(-1,5)
ax4.set_xlim(-1,5)
ax5.set_xlim(-1,5)

ax1.set_xlabel('deg/km', size = 12)
ax2.set_xlabel('deg/km',size = 12)
ax3.set_xlabel('deg/km',size = 12)
ax4.set_xlabel('deg/km',size = 12)
ax5.set_xlabel('deg/km',size = 12)

ax1.set_ylabel('Altitude (km)', size = 12)
ax2.set_ylabel('Altitude (km)',size = 12)
ax3.set_ylabel('Altitude (km)',size = 12)
ax4.set_ylabel('Altitude (km)',size = 12)
ax5.set_ylabel('Altitude (km)',size = 12)


ax1.set_title(' 8/26 00Z', size = 16)
ax2.set_title('8/26 12Z', size = 16)
ax3.set_title(' 8/27 00Z', size = 16)
ax4.set_title('8/27 12Z', size = 16)
ax5.set_title('8/28 00Z', size = 16)
subplots_adjust(wspace = 0.4, hspace = 0.3)

fig2 = gcf()
st = fig2.suptitle("Houston $K_{DP}$", size=24)
plt.subplots()


####Zh

fig2 = plt.figure(figsize = (18,12))
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)

ax1.plot(zh2,alt, color = 'b')
ax2.plot(zh3, alt, color = 'b')
ax3.plot(zh4, alt, color = 'b')
ax4.plot(zh5, alt, color = 'b')
ax5.plot(zh6, alt, color = 'b')

ax1.set_ylim(1,10)
ax2.set_ylim(1,10)
ax3.set_ylim(1,10)
ax4.set_ylim(1,10)
ax5.set_ylim(1,10)

ax1.set_xlabel('dBZ', size = 12)
ax2.set_xlabel('dBZ',size = 12)
ax3.set_xlabel('dBZ',size = 12)
ax4.set_xlabel('dBZ',size = 12)
ax5.set_xlabel('dBZ',size = 12)

ax1.set_ylabel('Altitude (km)', size = 12)
ax2.set_ylabel('Altitude (km)',size = 12)
ax3.set_ylabel('Altitude (km)',size = 12)
ax4.set_ylabel('Altitude (km)',size = 12)
ax5.set_ylabel('Altitude (km)',size = 12)


ax1.set_title(' 8/26 00Z', size = 16)
ax2.set_title('8/26 12Z', size = 16)
ax3.set_title(' 8/27 00Z', size = 16)
ax4.set_title('8/27 12Z', size = 16)
ax5.set_title('8/28 00Z', size = 16)
subplots_adjust(wspace = 0.4, hspace = 0.3)

fig2 = gcf()
st = fig2.suptitle("Houston $Z_{H}$", size=24)
plt.subplots()