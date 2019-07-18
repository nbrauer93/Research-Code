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
#from mpl_toolkits.basemap import Basemap,maskoceans,interp,shiftgrid
from scipy import io
from scipy import stats,signal
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator 
from matplotlib.colors import LinearSegmentedColormap
import glob
#import pyart
from scipy.interpolate import interp2d
from datetime import datetime
from matplotlib.dates import DateFormatter



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
alt = nc.variables['Altitude'][:]




    

files2 = glob.glob('nexrad*.nc')[:]    

zdr = np.ones((1441,28,4,4))*np.nan
kdp = np.ones((1441,28,4,4))*np.nan
zh = np.ones((1441,28,4,4))*np.nan
cc = np.ones((1441,28,4,4))*np.nan

zdrcrp = np.ones((1441,28,4,4))*np.nan
kdpcrp = np.ones((1441,28,4,4))*np.nan
zhcrp = np.ones((1441,28,4,4))*np.nan
cccrp = np.ones((1441,28,4,4))*np.nan


zdrlib = np.ones((1441,28,4,4))*np.nan
kdplib = np.ones((1441,28,4,4))*np.nan
zhlib = np.ones((1441,28,4,4))*np.nan
cclib = np.ones((1441,28,4,4))*np.nan



##Take fields at a single point (Houston) (227,365)

#Corpus Christi (27.8021,-97.4062) (134,268)


for n,j in enumerate(files2):
    data = read_file(j)
    print(j)
    zdr[n,:,:,:] = data['zdr']['values'][:,225:229,363:367]
    zh[n,:,:,:] = data['Z_H']['values'][:,225:229,363:367]
    cc[n,:,:,:] = data['cc']['values'][:,225:229,363:367]
    kdp[n,:,:,:] = data['kdp']['values'][:,225:229,363:367]

    
    #Liberty, TX: (30.0554,-94.7932) [:,242,393]
    zdrlib[n,:,:,:] = data['zdr']['values'][:,240:244,391:395]
    zhlib[n,:,:,:] = data['Z_H']['values'][:,240:244,391:395]
    cclib[n,:,:,:] = data['cc']['values'][:,240:244,391:395]
    kdplib[n,:,:,:] = data['kdp']['values'][:,240:244,391:395]
    
#%%

###Average over the five gridpoints by lat-lon


T,Z,I,J = zdr.shape
zdr_squish = zdr.reshape(T,Z,I*J, order = 'F')
kdp_squish = kdp.reshape(T,Z,I*J, order = 'F')
zh_squish = zh.reshape(T,Z,I*J, order = 'F')
cc_squish = cc.reshape(T,Z,I*J, order = 'F')



zhlib_squish = zhlib.reshape(T,Z,I*J, order = 'F')
zdrlib_squish = zdrlib.reshape(T,Z,I*J, order = 'F')
kdplib_squish = kdplib.reshape(T,Z,I*J, order = 'F')
cclib_squish = cclib.reshape(T,Z,I*J, order = 'F')




    
zdr_avg = np.ones((1441,28))*np.nan
kdp_avg = np.ones((1441,28))*np.nan    
zh_avg = np.ones((1441,28))*np.nan    
cc_avg = np.ones((1441,28))*np.nan 


zdr_avg_lib = np.ones((1441,28))*np.nan
kdp_avg_lib = np.ones((1441,28))*np.nan
zh_avg_lib = np.ones((1441,28))*np.nan
cc_avg_lib = np.ones((1441,28))*np.nan




for i in range(zdr.shape[0]):
    zdr_avg[i,:] = np.nanmean(zdr_squish[i,:,:],axis = 1) 
    zh_avg[i,:] = np.nanmean(zh_squish[i,:,:],axis = 1) 
    kdp_avg[i,:] = np.nanmean(kdp_squish[i,:,:],axis = 1) 
    cc_avg[i,:] = np.nanmean(cc_squish[i,:,:],axis = 1) 
    
    zdr_avg_lib[i,:] = np.nanmean(zdrlib_squish[i,:,:],axis = 1) 
    zh_avg_lib[i,:] = np.nanmean(zhlib_squish[i,:,:],axis = 1) 
    kdp_avg_lib[i,:] = np.nanmean(kdplib_squish[i,:,:],axis = 1) 
    cc_avg_lib[i,:] = np.nanmean(cclib_squish[i,:,:],axis = 1) 
 
    
#%%    
    
#Filter out scatterers associated with reflectivity = 0
    
'''    
no_scatterer = np.where((zh_avg>0))[0]    

zdr_avg2 = zdr_avg[no_scatterer]
kdp_avg2 = kdp_avg[no_scatterer]

zdr_avg_lib2 = zdr_avg_lib[no_scatterer]
kdp_avg_lib2 = kdp_avg_lib[no_scatterer]
'''

'''
zdr_avg_nan = np.ones((1441,28))*np.nan
kdp_avg_nan = np.ones((1441,28))*np.nan
zdr_avg_nan_lib = np.ones((1441,28))*np.nan
kdp_avg_nan_lib = np.ones((1441,28))*np.nan


for i in range(zdr_avg_nan.shape[0]):
    for j in range(zdr_avg_nan.shape[1]):
        if zh_avg[i,j] == 0:
            zdr_avg_nan[i,j] = np.nan
            kdp_avg_nan[i,j] = np.nan
            
        else:
            zdr_avg_nan[i,j] = zdr_avg[i,j]
            kdp_avg_nan[i,j] = kdp_avg[i,j]
            
            
        if zh_avg_lib[i,j] == 0:
             zdr_avg_nan[i,j] = np.nan
             kdp_avg_nan[i,j] = np.nan
        
        else:
            zdr_avg_nan_lib[i,j] = zdr_avg_lib[i,j]
            kdp_avg_nan_lib[i,j] = kdp_avg_lib[i,j]
            
'''            
            
             
    
            


zdr_avg_nan = np.ones((1441,28))*np.nan
kdp_avg_nan = np.ones((1441,28))*np.nan
zdr_avg_nan_lib = np.ones((1441,28))*np.nan
kdp_avg_nan_lib = np.ones((1441,28))*np.nan


for i in range(zdr_avg_nan.shape[0]):
    for j in range(zdr_avg_nan.shape[1]):
        if zdr_avg[i,j] == 0:
            zdr_avg_nan[i,j] = np.nan
            
            
        else:
            zdr_avg_nan[i,j] = zdr_avg[i,j]
            
            
            
        if zdr_avg_lib[i,j] == 0:
             zdr_avg_nan_lib[i,j] = np.nan
            
        
        else:
            zdr_avg_nan_lib[i,j] = zdr_avg_lib[i,j]
            
        
        if kdp_avg[i,j] == 0:
            kdp_avg_nan[i,j] = np.nan
            
        else:
            kdp_avg_nan[i,j] = kdp_avg[i,j]
            
            
        if kdp_avg_lib[i,j] == 0:
            kdp_avg_nan_lib[i,j] = np.nan
            
        else:
            kdp_avg_nan_lib[i,j] = kdp_avg_lib[i,j]
            
            











    
    
#%%   
 
#Time is a list of all times    
     
time = []    


for j in (files2):
    data2 = read_file(j)
    print(j)
    timejunk = data2['Analysis_time']
    time.append(timejunk)
    
    
#time2 = time[::3]
#time2 = time[::12]
#timelist = np.arange(0,480,1)
#timelist = np.arange(0,120,1)




altlist = np.array([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.])





    
    
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


zhcolor = reflect_ncdc()


from matplotlib.colors import ListedColormap
clevs = np.linspace(0.,2.0,20)
#colors = ['firebrick','red','darkorange','gold','green','lightgreen','lawngreen','aquamarine', 'cyan', 'deepskyblue','dodgerblue','blue', 'magenta']
colors = ['white','red','darkorange','gold','green','lightgreen','lawngreen','aquamarine', 'cyan', 'deepskyblue','dodgerblue','blue', 'magenta']
cmap_kdp = ListedColormap(colors)



clevs = np.linspace(0,4,20)
colors_zdr = ['lightgray','mediumorchid','darkblue','b','deepskyblue','aqua','springgreen','lime','chartreuse','greenyellow','yellow','orange','orangered','red']
cmap_zdr = ListedColormap(colors_zdr)


##Compute 20-minute running mean of radar data at each level

'''
zdravg = np.ones((480,28))*np.nan
zhavg = np.ones((480,28))*np.nan
kdpavg = np.ones((480,28))*np.nan
ccavg =  np.ones((480,28))*np.nan


zdravglib = np.ones((480,28))*np.nan
zhavglib = np.ones((480,28))*np.nan
kdpavglib = np.ones((480,28))*np.nan


zdravgcrp = np.ones((480,28))*np.nan
zhavgcrp = np.ones((480,28))*np.nan
kdpavgcrp = np.ones((480,28))*np.nan
ccavgcrp =  np.ones((480,28))*np.nan

for i in range(0,480):
    zdravg[i,:] = np.nanmean(zdr_avg[3*i:(3*i)+3,:], axis = 0)
    zhavg[i,:] = np.nanmean(zh_avg[3*i:(3*i)+3,:], axis = 0)
    kdpavg[i,:] = np.nanmean(kdp_avg[3*i:(3*i)+3,:], axis = 0)
    #ccavg[i,:] = np.nanmean(ccavg[3*i:(3*i)+3,:], axis = 0)
    
    zdravglib[i,:] = np.nanmean(zdr_avg_lib[3*i:(3*i)+3,:], axis = 0)
    zhavglib[i,:] = np.nanmean(zh_avg_lib[3*i:(3*i)+3,:], axis = 0)
    kdpavglib[i,:] = np.nanmean(kdp_avg_lib[3*i:(3*i)+3,:], axis = 0)
'''
   


    

    
    
#%%    
#HOUSTON 

time_cut = [times[8:] for times in time]

###Do this for reflectivity    
time2627 = time[0:288]
time_cut2627 = [times[11:] for times in time2627] 
time2728 = time[288:576]
time_cut2728 = [times[11:] for times in time2728] 
time2829 = time[576:864]
time_cut2829 = [times[11:] for times in time2829] 
time2930 = time[864:1152]
time_cut2930 = [times[11:] for times in time2930] 



'''
fig = plt.figure(figsize = [80,15])  
plt.ylabel('Altitude (km)', size = 14)
plt.ylim(1,15) 
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{H}$ (dBZ)', size = 20)  
c=plt.imshow(zhavg.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
fig.colorbar(c, orientation='horizontal', label = 'dBZ', fraction=0.046)
plt.xticks(np.arange(0, 480, step=80),time2[::80], rotation = 45)
plt.yticks(altlist[::3],altlist[::3])
'''
zh2627 = zh[0:288,:]
zh2728 = zh[288:576,:]
zh2829 = zh[576:864,:]

zh2627new = zh_avg[0:288,:]
zh2728new = zh_avg[288:576,:]
zh2829new = zh_avg[576:864,:]


alt2627new = alt[0:288]
alt2728new = alt[288:576]
alt2829new = alt[576:864]



'''    

fig = plt.figure(figsize = [20,5]) 
#plt.rcParams["figure.figsize"] = [10,20]
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{H}$  8/26-8/27', size = 28)  
k=plt.imshow(zh2627new.T,cmap=zhcolor,  interpolation = 'hanning',vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
fig.colorbar(k, orientation='horizontal', label = 'dBZ', fraction=0.1)
#plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
#plt.yticks(altlist[::6],altlist[::6])
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72])

'''

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
plt.title(r'Houston $Z_{H}$  8/26-8/27', size = 28) 
k = plt.contourf(zh2627new.T, levels = alt, cmap = zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
#plt.gca().invert_yaxis()
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72])
plt.show()






fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{H}$  8/27-8/28', size = 28)  
k=plt.imshow(zh2728new.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dBZ', fraction=0.03)
plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{H}$  8/28-8/29', size = 28)  
k=plt.imshow(zh2829new.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
#plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2829[::72])






#And now for differential reflectivity
'''
fig = plt.figure(figsize = [80,15])  
plt.ylabel('Altitude (km)', size = 14)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{DR}$ ', size = 20)  
h=plt.imshow(zdravg.T, cmap=zdrcolor, interpolation = 'hanning')
plt.gca().invert_yaxis()
fig.colorbar(h, orientation='horizontal', label = 'dB',fraction=0.046)
plt.xticks(np.arange(0, 480, step=80),time2[::80], rotation = 45)
plt.yticks(altlist[::3],altlist[::3])
'''

zdr2627new = zdr_avg_nan[0:288,:]
zdr2728new = zdr_avg[288:576,:]
zdr2829new = zdr_avg[576:864,:]



zdr2627 = zdr[0:288,:]
zdr2728 = zdr[288:576,:]
zdr2829 = zdr[576:864,:]



fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{DR}$  8/26-8/27', size = 28)  
k=plt.imshow(zdr2627new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72])



fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{DR}$  8/27-8/28', size = 28)  
k=plt.imshow(zdr2728new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{DR}$  8/28-8/29', size = 28)  
k=plt.imshow(zdr2829new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2829[::72])



#Specific Differential Phase
'''
fig = plt.figure(figsize = [80,15])  
plt.ylabel('Altitude (km)', size = 14)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $K_{DP}$ (deg/km)', size = 20)  
k=plt.imshow(kdpavg.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.046)
plt.xticks(np.arange(0, 480, step=80),time2[::80], rotation = 45)
plt.yticks(altlist[::3],altlist[::3])
'''

kdp2627new = kdp_avg[0:288,:]
kdp2728new = kdp_avg[288:576,:]
kdp2829new = kdp_avg[576:864,:]

kdp2627 = kdp[0:288,:]
kdp2728 = kdp[288:576,:]
kdp2829 = kdp[576:864,:]

#Chop it into times 26-27, 27-28, 28-29
fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $K_{DP}$  8/26-8/27', size = 28)  
k=plt.imshow(kdp2627new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72], rotation = 45)

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $K_{DP}$  8/27-8/28', size = 28)  
k=plt.imshow(kdp2728new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $K_{DP}$  8/28-8/29', size = 28)  
k=plt.imshow(kdp2829new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])



#Correlation coefficient
'''
fig = plt.figure(figsize = [80,15])  
plt.ylabel('Altitude (km)', size = 14)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $\rho_{hv}$', size = 20)  
l=plt.imshow(ccavg.T, cmap= zhcolor , interpolation = 'hanning', vmin = 0.90, vmax = 1.01)
plt.gca().invert_yaxis()
fig.colorbar(l, orientation='horizontal',fraction=0.046)
plt.xticks(np.arange(0, 480, step=80),time2[::80])
plt.yticks(altlist[::3],altlist[::3])

'''




#%%

#########For Liberty, TX

zhlib2627 = zhlib[0:288,:]
zhlib2728 = zhlib[288:576,:]
zhlib2829 = zhlib[576:864,:]
zhlib2930 = zhlib[864:1152,:]

zhlib2627new = zh_avg_lib[0:288,:]
zhlib2728new = zh_avg_lib[288:576,:]
zhlib2829new = zh_avg_lib[576:864,:]
zhlib2930new = zh_avg_lib[864:1152,:]


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{H}$  8/26-8/27', size = 28)  
k=plt.imshow(zhlib2627new.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dBZ', fraction=0.03)
plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72])



fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{H}$  8/27-8/28', size = 28)  
k=plt.imshow(zhlib2728new.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dBZ', fraction=0.03)
plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{H}$  8/28-8/29', size = 28)  
k=plt.imshow(zhlib2829new.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dBZ', fraction=0.03)
plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2829[::72])


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{H}$  8/29-8/30', size = 28)  
k=plt.imshow(zhlib2930new.T, cmap=zhcolor, interpolation = 'hanning', vmin = 0, vmax = 75)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dBZ', fraction=0.03)
plt.colorbar().set_label(label='dBZ',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72])




#And now for differential reflectivity
'''
fig = plt.figure(figsize = [80,15])  
plt.ylabel('Altitude (km)', size = 14)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $Z_{DR}$ ', size = 20)  
h=plt.imshow(zdravg.T, cmap=zdrcolor, interpolation = 'hanning')
plt.gca().invert_yaxis()
fig.colorbar(h, orientation='horizontal', label = 'dB',fraction=0.046)
plt.xticks(np.arange(0, 480, step=80),time2[::80], rotation = 45)
plt.yticks(altlist[::3],altlist[::3])
'''



zdrlib2627 = zdrlib[0:288,:]
zdrlib2728 = zdrlib[288:576,:]
zdrlib2829 = zdrlib[576:864,:]
zdrlib2930 = zdrlib[864:1152,:]

zdrlib2627new = zdr_avg_lib[0:288,:]
zdrlib2728new = zdr_avg_lib[288:576,:]
zdrlib2829new = zdr_avg_lib[576:864,:]
zdrlib2930new = zdr_avg_lib[864:1152,:]

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{DR}$  8/26-8/27', size = 28)  
k=plt.imshow(zdrlib2627new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72])



fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{DR}$  8/27-8/28', size = 28)  
k=plt.imshow(zdrlib2728new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{DR}$  8/28-8/29', size = 28)  
k=plt.imshow(zdrlib2829new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2829[::72])

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $Z_{DR}$  8/29-8/30', size = 28)  
k=plt.imshow(zdrlib2930new.T, cmap=zdrcolor, interpolation = 'hanning', vmin = -3, vmax = 4)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'dB', fraction=0.03)
plt.colorbar().set_label(label='dB',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2930[::72])


#Specific Differential Phase
'''
fig = plt.figure(figsize = [80,15])  
plt.ylabel('Altitude (km)', size = 14)
#plt.xlabel('Time', size = 14)
plt.title(r'Houston $K_{DP}$ (deg/km)', size = 20)  
k=plt.imshow(kdpavg.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.046)
plt.xticks(np.arange(0, 480, step=80),time2[::80], rotation = 45)
plt.yticks(altlist[::3],altlist[::3])
'''

kdplib2627 = kdplib[0:288,:]
kdplib2728 = kdplib[288:576,:]
kdplib2829 = kdplib[576:864,:]
kdplib2930 = kdplib[864:1152,:]

kdplib2627new = kdp_avg_lib[0:288,:]
kdplib2728new = kdp_avg_lib[288:576,:]
kdplib2829new = kdp_avg_lib[576:864,:]
kdplib2930new = kdp_avg_lib[864:1152,:]



#Chop it into times 26-27, 27-28, 28-29
fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $K_{DP}$  8/26-8/27', size = 28)  
k=plt.imshow(kdplib2627new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2627[::72], rotation = 45)

fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $K_{DP}$  8/27-8/28', size = 28)  
k=plt.imshow(kdplib2728new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $K_{DP}$  8/28-8/29', size = 28)  
k=plt.imshow(kdplib2829new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])


fig = plt.figure(figsize = [20,5]) 
plt.ylabel('Altitude (km)', size = 16)
#plt.xlabel('Time', size = 14)
plt.title(r'Liberty $K_{DP}$  8/29-8/30', size = 28)  
k=plt.imshow(kdplib2930new.T, cmap=cmap_kdp, interpolation = 'hanning', vmin = 0, vmax = 2)
plt.gca().invert_yaxis()
#fig.colorbar(k, orientation='horizontal', label = 'deg/km', fraction=0.03)
plt.colorbar().set_label(label='deg/km',size=18,weight='bold')
plt.yticks(altlist[::6],altlist[::6])
plt.xticks(np.arange(0,288, step = 72),time_cut2728[::72])


    
    

              
              