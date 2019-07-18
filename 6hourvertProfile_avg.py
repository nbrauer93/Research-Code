#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:41:48 2019

@author: noahbrauer
"""

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
    
    
    



file = 'nexrad_3d_v4_0_20170828T211500Z.nc'

nc = Dataset(file ,'r')

lat = nc.variables['Latitude'][:]
lon = nc.variables['Longitude'][:] - 360
alt = nc.variables['Altitude'][:]





files2 = glob.glob('nexrad*.nc')[:]

zdr = np.ones((1441,28,4,4))*np.nan
kdp = np.ones((1441,28,4,4))*np.nan
cc = np.ones((1441,28,4,4))*np.nan
zh = np.ones((1441,28,4,4))*np.nan

###3 km radar fields

for n,j in enumerate(files2):
    print(j)
    data = read_file(j)

    zdr[n,:] = data['zdr']['values'][:,225:229,363:367]
    kdp[n,:] = data['kdp']['values'][:,225:229,363:367]
    zh[n,:] = data['Z_H']['values'][:,225:229,363:367]
    cc[n,:] = data['cc']['values'][:,225:229,363:367]  

#%%

T,Z,I,J = zdr.shape
zdr_squish = zdr.reshape(T,Z,I*J, order = 'F')
zh_squish = zh.reshape(T,Z,I*J, order = 'F')
kdp_squish = kdp.reshape(T,Z,I*J, order = 'F')
cc_squish = cc.reshape(T,Z,I*J, order = 'F')


#Filter out values of zh,zdr,kdp,rhohv where zh = 0

zh_nan = np.ones((1441,28,16))*np.nan
zdr_nan = np.ones((1441,28,16))*np.nan
kdp_nan = np.ones((1441,28,16))*np.nan
cc_nan = np.ones((1441,28,16))*np.nan


for i in range(zh_nan.shape[0]):
    for j in range(zh_nan.shape[1]):
        for k in range(zh_nan.shape[2]):
            
            if cc_squish[i,j,k] <0.8:
                zh_nan[i,j,k] = np.nan
                zdr_nan[i,j,k] = np.nan
                kdp_nan[i,j,k] = np.nan
                cc_nan[i,j,k] = np.nan
                
            else:
                zh_nan[i,j,k] = zh_squish[i,j,k]
                zdr_nan[i,j,k] = zdr_squish[i,j,k]
                kdp_nan[i,j,k] = kdp_squish[i,j,k]
                cc_nan[i,j,k] = cc_squish[i,j,k]
                



zdravg = np.ones((1441,28))*np.nan
zhavg = np.ones((1441,28))*np.nan
kdpavg = np.ones((1441,28))*np.nan
ccavg = np.ones((1441,28))*np.nan

#Take spatial mean for each time step

for i in range(zh.shape[0]):
    zdravg[i,:] = np.nanmean(zdr_squish[i,:,:], axis = 1)
    zhavg[i,:] = np.nanmean(zh_squish[i,:,:], axis = 1)
    kdpavg[i,:] = np.nanmean(kdp_squish[i,:,:], axis = 1)
    ccavg[i,:] = np.nanmean(cc_squish[i,:,:], axis = 1)
    
    
#Now take 6-hour means
    
    

vertical_zh = np.ones((20,28))*np.nan    
vertical_zdr = np.ones((20,28))*np.nan   
vertical_kdp = np.ones((20,28))*np.nan   
vertical_cc = np.ones((20,28))*np.nan   


for i in range(vertical_zh.shape[0]):
    for j in range(vertical_zh.shape[1]):
        vertical_zh[i,j] = np.nanmean(zhavg[72*i:72*i+72,j], axis = 0)
        vertical_zdr[i,j] = np.nanmean(zdravg[72*i:72*i+72,j], axis = 0)
        vertical_kdp[i,j] = np.nanmean(kdpavg[72*i:72*i+72,j], axis = 0)
        vertical_cc[i,j] = np.nanmean(ccavg[72*i:72*i+72,j], axis = 0)


#%%
#Now plot
altlist = np.array([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.])
fig = plt.figure(figsize = [10,6])
plt.subplots_adjust(bottom=0.1, right=1.5, top=2)
    

plt.subplot(2, 2, 1)
plt.xlabel('dBZ', size = 18)
plt.ylabel('Altitude (km)', size = 18)
plt.xlim(0,35)
plt.ylim(np.min(altlist),np.max(altlist))
plt.yticks(altlist[::3],altlist[::3])
plt.title(r'Houston 6 hour mean $Z_{H}$ 8/28-8/29', size = 25)
plt.plot(vertical_zh[8,:], alt, color = 'r', label = '00Z-06Z', linestyle = 'dashed')
plt.plot(vertical_zh[9,:], alt, color = 'g', label = '06Z-12Z')
plt.plot(vertical_zh[10,:], alt, color = 'b', label = '12Z-18Z')
plt.plot(vertical_zh[11,:], alt, color = 'k', label = '18Z-00Z', linestyle = 'dashed')
plt.legend()



plt.subplot(2, 2, 2)
plt.xlabel('dB', size = 18)
plt.ylabel('Altitude (km)', size = 18)
plt.xlim(0,1.5)
plt.ylim(np.min(altlist),np.max(altlist))
plt.yticks(altlist[::3],altlist[::3])
plt.title(r'Houston 6 hour mean $Z_{DR}$ 8/28-8/29', size = 25)
plt.plot(vertical_zdr[8,:], alt, color = 'r', label = '00Z-06Z', linestyle = 'dashed')
plt.plot(vertical_zdr[9,:], alt, color = 'g', label = '06Z-12Z')
plt.plot(vertical_zdr[10,:], alt, color = 'b', label = '12Z-18Z')
plt.plot(vertical_zdr[11,:], alt, color = 'k', label = '18Z-00Z', linestyle = 'dashed')
plt.legend()

plt.subplot(2, 2, 3)
plt.xlabel('deg/km', size = 18)
plt.ylabel('Altitude (km)', size = 18)
plt.xlim(-0.25,1)
plt.ylim(np.min(altlist),np.max(altlist))
plt.yticks(altlist[::3],altlist[::3])
plt.title(r'Houston 6 hour mean $K_{DP}$ 8/28-8/29', size = 25)
plt.plot(vertical_kdp[8,:], alt, color = 'r', label = '00Z-06Z', linestyle = 'dashed')
plt.plot(vertical_kdp[9,:], alt, color = 'g', label = '06Z-12Z')
plt.plot(vertical_kdp[10,:], alt, color = 'b', label = '12Z-18Z')
plt.plot(vertical_kdp[11,:], alt, color = 'k', label = '18Z-00Z', linestyle = 'dashed')
plt.legend()

plt.subplot(2, 2, 4)
plt.xlabel(r'$\rho_{hv}$',size = 18)
plt.ylabel('Altitude (km)', size = 18)
plt.xlim(0,1)
plt.ylim(np.min(altlist),np.max(altlist))
plt.yticks(altlist[::3],altlist[::3])
plt.title(r'Houston 6 hour mean $\rho_{hv}$ 8/28-8/29', size = 25)
plt.plot(vertical_cc[8,:], alt, color = 'r', label = '00Z-06Z', linestyle = 'dashed')
plt.plot(vertical_cc[9,:], alt, color = 'g', label = '06Z-12Z')
plt.plot(vertical_cc[10,:], alt, color = 'b', label = '12Z-18Z')
plt.plot(vertical_cc[11,:], alt, color = 'k', label = '18Z-00Z', linestyle = 'dashed')
plt.legend()
  

plt.show()    






    
