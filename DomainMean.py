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





files2 = glob.glob('nexrad*.nc')[:]

zdr = np.ones((1440,528,672))*np.nan
kdp = np.ones((1440,528,672))*np.nan
cc = np.ones((1440,528,672))*np.nan
zh = np.ones((1440,528,672))*np.nan

###3 km radar fields

for n,j in enumerate(files2):
    print(j)
    data = read_file(j)

    zdr[n,:,:] = data['zdr']['values'][4,:,:]
    kdp[n,:,:] = data['kdp']['values'][4,:,:]
    zh[n,:,:] = data['Z_H']['values'][4,:,:]
    cc[n,:,:] = data['cc']['values'][4,:,:]    
    
    
    
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
##Houston Domain

xlim = np.array([28.,31.]); ylim = np.array([-97.,-94.])

ilat = np.where((lat>=np.min(xlim))&(lat<=np.max(xlim)))[0]
ilon = np.where((lon>=np.min(ylim))&(lon<=np.max(ylim)))[0]

'''
ilat = np.where((lat>=xlim[0])& (lat<=xlim[1]))[0]
ilon = np.where((lon>=ylim[0])& (lon<=ylim[1]))[0]
'''
#Define a grid using these points

lat2, lon2 = np.meshgrid(lat[ilat], lon[ilon])

zdrhou = zdr[:,ilat[0]:ilat[-1]+1,ilon[0]:ilon[-1]+1]
cchou = cc[:,ilat[0]:ilat[-1]+1,ilon[0]:ilon[-1]+1]
zhhou = zh[:,ilat[0]:ilat[-1]+1,ilon[0]:ilon[-1]+1]
kdphou = kdp[:,ilat[0]:ilat[-1]+1,ilon[0]:ilon[-1]+1]

#Take mean of each polarimetric radar field over the time axis

T,I,J = zdrhou.shape
zdrsquish = zdrhou.reshape(T,I*J)
zhsquish = zhhou.reshape(T,I*J)
kdpsquish = kdphou.reshape(T,I*J)



zdrmeanint = np.ones((9,20736))*np.nan
kdpmeanint = np.ones((9,20736))*np.nan
zhmeanint = np.ones((9,20736))*np.nan

for i in range(8):
    for j in range(20736):
        zdrmeanint[i,j] = np.nanmean(zdrsquish[576+72*(i-1):576+72*i,j], axis = 0)
        kdpmeanint[i,j] = np.nanmean(kdpsquish[576+72*(i-1):576+72*i,j], axis = 0)
        zhmeanint[i,j] = np.nanmean(zhsquish[576+72*(i-1):576+72*i,j], axis = 0)







zdrmean2 = zdrmeanint.reshape(9,144,144)
kdpmean2 = kdpmeanint.reshape(9,144,144)
zhmean2 = zhmeanint.reshape(9,144,144)








xlim2 = np.array([25.,28.]); ylim2 = np.array([-99.,-96.])

ilat2 = np.where((lat>=np.min(xlim2))&(lat<=np.max(xlim2)))[0]
ilon2 = np.where((lon>=np.min(ylim2))&(lon<=np.max(ylim2)))[0]




lat3, lon3 = np.meshgrid(lat[ilat2], lon[ilon2]) #For Corpus


    #Corpus Christi Domain
zdrcrp = zdr[:,ilat2[0]:ilat2[-1]+1,ilon2[0]:ilon2[-1]+1]
cccrp = cc[:,ilat2[0]:ilat2[-1]+1,ilon2[0]:ilon2[-1]+1]
zhcrp = zh[:,ilat2[0]:ilat2[-1]+1,ilon2[0]:ilon2[-1]+1]
kdpcrp = kdp[:,ilat2[0]:ilat2[-1]+1,ilon2[0]:ilon2[-1]+1]

    #Corpus Christi
zdrsquish2 = zdrcrp.reshape(T,I*J)
zhsquish2 = zhcrp.reshape(T,I*J)
kdpsquish2 = kdpcrp.reshape(T,I*J)


zdrmeanintcrp = np.ones((9,20736))*np.nan
kdpmeanintcrp = np.ones((9,20736))*np.nan
zhmeanintcrp = np.ones((9,20736))*np.nan

for i in range(0,8):
    for j in range(20736):
        zdrmeanintcrp[i,j] = np.nanmean(zdrsquish2[72*i:72*i+72,j], axis = 0)
        kdpmeanintcrp[i,j] = np.nanmean(kdpsquish2[72*i:72*i+72,j], axis = 0)
        zhmeanintcrp[i,j] = np.nanmean(zhsquish2[72*i:72*i+72,j], axis = 0)


    #Corpus Christi averages
zdrmeancrp = zdrmeanintcrp.reshape(9,144,144)
kdpmeancrp = kdpmeanintcrp.reshape(9,144,144)
zhmeancrp = zhmeanintcrp.reshape(9,144,144)

#%%


cmin = -1.; cmax = 3.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(18,12))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-97,-94]); ylim = np.array([28,31])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='h')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon2,lat2,zdrmean2[7,:,:],clevs,cmap=zdrcolor,extend='both') #plot lat, lon, and North Pacific SST Anomalies

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

plt.title(r'Mean 3km $Z_{DR}$ 8/27 06Z-12Z',name='Calibri',size=16)
plt.show(block=False)

#%%


cmin = 0.; cmax = 80.; cint = 5.; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    #plt.figure()
plt.figure(figsize=(18,12))
    #cornpmm, lon = shiftgrid(180., corrnpmm, lon, start=False)
    #lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, 90, 180))
xlim = np.array([-97,-94]); ylim = np.array([28,31])
#xlim = np.array([-110,-85]); ylim = np.array([10,50])
    #parallels = np.arange(27.,35.,1.)
    # labels = [left,right,top,bottom]
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='h')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  #m.drawparallels(parallels,labels=[True,True,True,True], size = 20) #m.fillcontinents(color='Black');m.drawparallels(np.arange(-180.,180.,30.), labels = [1,0,0,0]);m.drawmeridians(np.arange(-180,180,30),labels = [0,0,0,1])
    #xCoord,yCoord = m(lat2, lon2)  ###Add lat-lon here
cs = m.contourf(lon2,lat2,zhmean2[7,:,:],clevs,cmap=zhcolor,extend='both') #plot lat, lon, and North Pacific SST Anomalies

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

plt.title(r'Mean 3km $Z_{H}$ 8/27 06Z-12Z',name='Calibri',size=16)
plt.show(block=False)

