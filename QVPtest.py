#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:15:56 2018

@author: noahbrauer
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 19:04:00 2018

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
from matplotlib.colors import LinearSegmentedColormap
import glob
import pyart
from scipy.interpolate import interp2d


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
    wcc = id.variables['wCorrelationCoefficient']
    
    
    
     

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
           'cc'                      : cc,  \
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
nc = Dataset(file, 'r')

lon = nc.variables['Longitude'][:] -360
lat = nc.variables['Latitude'][:] 
alt = nc.variables['Altitude'][0:15] #1-10 km


    

files2 = glob.glob('nexrad*.nc')[:]    

zdr = np.ones((1441,15))*np.nan
kdp = np.ones((1441,15))*np.nan
zh = np.ones((1441,15))*np.nan
cc = np.ones((1441,15))*np.nan


##Take fields at a single point (Houston) from 1 - 10 km AGL


for n,j in enumerate(files2):
    data = read_file(j)
    print(j)
    zdr[n,:] = data['zdr']['values'][0:15,227,365]
    zh[n,:] = data['Z_H']['values'][0:15,227,365]
    cc[n,:] = data['cc']['values'][0:15,227,365]
    kdp[n,:] = data['kdp']['values'][0:15,227,365]
 
    
    
    
#%%

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




zhtime = zh[200:400,:]
zhtime2 = np.flip(zhtime,axis = 1)


'''
def Plot_qvp(radar_data,color):
    fig = plt.figure(figsize = [15,15])  
    plt.ylabel('Altitude (km)', size = 14)
    plt.xlabel('Time', size = 14)
    plt.title(r'$Z_{H}$ (dBZ)', size = 16)  
    plotted = plt.imshow(radar_data.T, cmap=color)
    
    return plotted
'''



    

#plt.ylim(0,10)    
fig = plt.figure(figsize = [15,15])  
plt.ylabel('Altitude (km)', size = 14)
plt.xlabel('Time', size = 14)
plt.title(r'$Z_{H}$ (dBZ)', size = 16)  
plt.imshow(zhtime2.T, cmap=zhcolor)



    
    
    
    
    
    

              
              