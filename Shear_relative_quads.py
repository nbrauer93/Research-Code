#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 13:45:03 2018

@author: noahbrauer
"""

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.basemap import Basemap,maskoceans,interp
from netCDF4 import Dataset, num2date, MFDataset
import pyart
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime





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

file = 'nexrad_3d_v4_0_20170826T000000Z.nc'
nc = Dataset(file, 'r')


radar_data = read_file(file)
zh = radar_data['Z_H']['values'][:]
radarlat = nc.variables['Latitude'][:]
radarlon = nc.variables['Longitude'][:]-360
alt = nc.variables['Altitude'][:]

radarlat2,radarlon2 = np.meshgrid(radarlat,radarlon)

track = 'Harvey_Best_Track.txt'

harvey = np.loadtxt(track, dtype = 'str', skiprows = 1)

#Lat-lon of actual hurricane location

latitude = harvey[:,4] 
longitude = harvey[:,5]
time = harvey[:,:2]
pres = harvey[:,7]


#Convert from degrees to km

def latlon_to_dist(lat_input2,lat_input1,lon_input2,lon_input1):
    distlat = (lat_input2-lat_input1)*111
    distlon = (lon_input2 - lon_input1)*111
    totaldist = np.sqrt((distlat)**2 + (distlon)**2)
    
    return(totaldist)
        
#These are the coordinates of Harvey
    
coords = np.array([latitude,longitude],dtype = 'float').T # 74x2
centerpoint = coords[41,:] #This is the centerpoint at the time of the radar image

################ Based off center, define eyewall, inner rainbands, and outer rainbands ###############

#First, take differences in lat-lon coordinates equating to 15,60, and 110 km 

sum_eyewall = 15/111  
sum_inner = 60/111
sum_outer = 110/111

#These are intervals from the centerpoint defining each annuli
eyewall_dist = np.where(coords[:,:]<=sum_eyewall+coords[:,:])[0]
inner_dist = np.where((coords[:,:]+sum_eyewall>coords[:,:])&(coords[:,:]>=sum_inner+coords[:,:]))[0]
outer_dist = np.where((coords[:,:]+sum_inner>coords[:,:])&(coords[:,:]>=sum_outer+coords[:,:]))[0]
    

#These variables 


eyewall_lat_top = np.ones((74,2))*np.nan
inner_lat_top = np.ones((74,2))*np.nan
outer_lat_top = np.ones((74,2))*np.nan
eyewall_lat_bot = np.ones((74,2))*np.nan
inner_lat_bot = np.ones((74,2))*np.nan
outer_lat_bot = np.ones((74,2))*np.nan



for i in range(coords.shape[0]):
    eyewall_lat_top[i,:] = coords[i,:]+sum_eyewall
    inner_lat_top[i,:] = coords[i,:]+sum_inner
    outer_lat_top[i,:] = coords[i,:] + sum_outer
    eyewall_lat_bot[i,:] = coords[i,:]-sum_eyewall
    inner_lat_bot[i,:] = coords[i,:]-sum_inner
    outer_lat_bot[i,:] = coords[i,:] - sum_outer


eyewall_lon_top = np.ones((74,2))*np.nan
inner_lon_top = np.ones((74,2))*np.nan
outer_lon_top = np.ones((74,2))*np.nan 
eyewall_lon_bot = np.ones((74,2))*np.nan
inner_lon_bot = np.ones((74,2))*np.nan
outer_lon_bot = np.ones((74,2))*np.nan 
    
for j in range(coords.shape[1]):
   eyewall_lon_top[:,j] = coords[:,j]+sum_eyewall
   inner_lon_top[:,j] = coords[:,j]+sum_inner
   outer_lon_top[:,j] = coords[:,j] + sum_outer   
   eyewall_lon_bot[:,j] = coords[:,j]-sum_eyewall
   inner_lon_bot[:,j] = coords[:,j]-sum_inner
   outer_lon_bot[:,j] = coords[:,j] - sum_outer    


#Now we have defined the lat-lon coordinates for each range from the center (coords); combine new coords into a single array
   

eyewall_lat_only_top = eyewall_lat_top[:,0]
inner_lat_only_top = inner_lat_top[:,0]
outer_lat_only_top = outer_lat_top[:,0]
eyewall_lon_only_top = eyewall_lon_top[:,1]
inner_lon_only_top = inner_lon_top[:,1]
outer_lon_only_top = outer_lon_top[:,1]

eyewall_lat_only_bot = eyewall_lat_bot[:,0]
inner_lat_only_bot = inner_lat_bot[:,0]
outer_lat_only_bot = outer_lat_bot[:,0]
eyewall_lon_only_bot = eyewall_lon_bot[:,1]
inner_lon_only_bot = inner_lon_bot[:,1]
outer_lon_only_bot = outer_lon_bot[:,1]

###Defining max and min points for lat and lon. All radar variables within these outer and inner ranges will be averaged over space. 

   
eyewall_annulus_max = np.array([eyewall_lat_only_top,eyewall_lon_only_top]).T
eyewall_annulus_min = np.array([eyewall_lat_only_bot,eyewall_lon_only_bot]).T
inner_annulus_max = np.array([inner_lat_only_top,inner_lon_only_top]).T
inner_annulus_min = np.array([inner_lat_only_bot,inner_lon_only_bot]).T
outer_annulus_max = np.array([outer_lat_only_top,outer_lon_only_top]).T
outer_annulus_min = np.array([outer_lat_only_bot,outer_lon_only_bot]).T


   

#Assign values of reflectivity that fall within each annulus

###Test for this case at 8/26 00Z (index 41)

###Extracts only the values that fall within the  assigned radius for each annulus
radar_eyewall_lat = np.where((eyewall_annulus_max[41,0]>=radarlat)&(eyewall_annulus_min[41,0]<=radarlat))[0]
radar_eyewall_lon = np.where((eyewall_annulus_max[41,1]>=radarlon)&(eyewall_annulus_min[41,1]<=radarlon))[0]
radar_inner_lat = np.where((inner_annulus_max[41,0]>=radarlat)&(inner_annulus_min[41,0]<=radarlat))[0]
radar_inner_lon = np.where((inner_annulus_max[41,1]>=radarlon)&(inner_annulus_min[41,1]<=radarlon))[0]
radar_outer_lat = np.where((outer_annulus_max[41,0]>=radarlat)&(outer_annulus_min[41,0]<=radarlat))[0]
radar_outer_lon = np.where((outer_annulus_max[41,1]>=radarlon)&(outer_annulus_min[41,1]<=radarlon))[0]




#Assign the indices to the radar data
zh_eyewall = zh[:,radar_eyewall_lat,radar_eyewall_lon]
zh_inner = zh[:,radar_inner_lat,radar_inner_lon]
zh_outer = zh[:,radar_outer_lat,radar_outer_lon]




#Take spatial mean of each annuli; will output a spatially-averaged vertical profile
zh_eyewall_mean = np.nanmean(zh_eyewall,axis=1)
zh_inner_mean = np.nanmean(zh_inner,axis=1)
zh_outer_mean = np.nanmean(zh_outer,axis=1)



##Set labels for altitude

alt_label = ['1','1.5','2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']



plt.figure(figsize=(10,6))
plt.plot(zh_eyewall_mean,alt_label, color = 'b', linestyle = 'dashed', label = 'Eyewall')
plt.plot(zh_inner_mean,alt_label, color = 'r', label = 'Inner Rainbands')
plt.plot(zh_outer_mean,alt_label, color = 'k', label = 'Outer Rainbands')
plt.xlabel('Reflectivity Factor (dBZ)',weight = 'bold', size = 16)
plt.ylabel('Altitude (km)', weight = 'bold', size = 16)
plt.title('Reflectivity Factor Spatial Mean 8/26 00Z', weight = 'bold', size = 20)
plt.legend()
plt.show()






lon_axis_coords = np.ones(160)*np.nan
lat_axis_coords = np.ones(160)*np.nan
'''
for i in range(np.min(outer_lon_top[:,0]),np.max(outer_lon_top[:,0])):
    for j in range(np.min(outer_lat_top[:,1]),np.max(outer_lat_top[:,1])):
        
'''        




plt.figure(figsize=(18,12))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-100,-92]); ylim = np.array([24,31])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(radarlon2,radarlat2,zh[4,:,:].T, cmap = zhcolor) #plot lat, lon, and North Pacific SST Anomalies
m.drawcounties()
###Plot center of Harvey from BT data on here: coord[47]
    
    
####These few lines plot the entire track of Harvey,not just a single point in time    
'''
xjunk,yjunk = m(coords[:,1],coords[:,0])
m.plot(xjunk,yjunk,'ko',markersize = 10)
plt.show()
'''
#################


'''
#x2star,y2star = m(centerpoint[1],centerpoint[0])
#m.plot(x2star,y2star,'ko',markersize=10)
x,y = m(centerpoint[0],centerpoint[1])
x2,y2 = m(centerpoint[0],centerpoint[0]+sum_inner)
c= create_circle((x,y),y2-y)
plot_circle(c)
plt.show() 
'''


#%%


ufile = 'uwnd.201708.nc'
vfile = 'vwnd.201708.nc'

uwind = Dataset(ufile, 'r')
vwind = Dataset(vfile, 'r')

narr = {}

time = uwind.variables['time'][:]
timeUnits = uwind.variables['time'].units
tmpDates = num2date(time,timeUnits,calendar='gregorian')
narr['date'] = np.asarray([datetime(d.year,d.month,d.day) for d in tmpDates])
narr['day'] = np.asarray([d.day for d in narr['date']])
narr['month'] = np.asarray([d.month for d in narr['date']])
narr['year'] = np.asarray([d.year for d in narr['date']])


aug_index = np.where((narr['day']>=25)&(narr['month']==8))[0]

u = uwind.variables['uwnd'][aug_index,:,:,:]
v = vwind.variables['vwnd'][aug_index,:,:,:]
level = uwind.variables['level'][:]

u850 = u[:,6,:,:]
v850 = v[:,6,:,:]
u200 = u[:,24,:,:]
v200 = v[:,24,:,:]

def wind_mag (u,v):
    wind = np.sqrt((u**2)+(v**2))
    return wind

wind850 = wind_mag(u850,v850)
wind200 = wind_mag(u200,v200)


#Compute the deep-layer (850-200 mb shear)

shear = np.abs(wind200-wind850)



###Okay cool. Now split into DL,DR,UL,UR qudrants. Will need the deep-layer shear vector here.

#Designate a center line along the center point;


