# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 11:10:52 2018

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
import pickle
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


for n,j in enumerate(files2):
    data = read_file(j)

    zdr[n,:,:] = data['zdr']['values'][4,::3,::3]
    kdp[n,:,:] = data['kdp']['values'][4,::3,::3]
    
    
    
    #Z, I, J = z_h.shape ##check that these are in the right order (z,y,x)
    
    
#    I2,J2 = zdr.shape
#    I3,J3 = kdp.shape
#    
#    #z_h2 = z_h.reshape((Z, I*J), order='F')
#    
#    zdr2 = zdr.reshape((I2*J2), order = 'F')
#    kdp2 = kdp.reshape((I3*J3), order = 'F')
#    
#    
#    #Z_H.append(z_h[0,:])
#    
#    Zdr.append(zdr2[:])
#    Kdp.append(kdp2[:])
    
    
#%%
    
ncFile = 'apcp.2017.nc'

narr = {}

with Dataset(ncFile,'r') as nc:
    lon = nc.variables['lon'][:]
    lon2 = lon+360
    lat = nc.variables['lat'][:]
    lam = nc.variables['Lambert_Conformal'][:]
    time = nc.variables['time'][:]
    precip = nc.variables['apcp'][:].T/25.4 #shape: 349/277/2920 = IJT convert from mm to in
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

precip3 = precip[:,nanremove] #Goes into EOF machine
precip_ano_no_nan = precip_ano[:,nanremove] #Goes into regression

    
#rprecip = np.ones((precip.shape[1],))*np.nan

#nancov(precip3)
    
    
    
    
    
    


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
a = np.concatenate((zdr3.reshape(T,I*J), kdp3.reshape(T,I*J),precip3), axis = 1)



#Import the lat-lon for radar

radarlon = data['x']['values'][::3]-360
radarlat = data['y']['values'][::3]

radarlat2,radarlon2 = np.meshgrid(radarlat,radarlon)


###Plot reflectivity showing location of Hurricane Harvey for reference









#%%

###Compute the EOFS


def nancov(A,B):
    A2 = A.copy(); B2 = B.copy()
    NmatA = np.where(np.isnan(A2),0,1).squeeze()
    NmatB = np.where(np.isnan(B2),0,1).squeeze()
    A2[np.isnan(A2)] = 0.; B2[np.isnan(B2)] = 0.
    C = np.dot(A2,B2)/np.dot(NmatA,NmatB)
    return C




def Compute_EOFs(anom,NMODES=None,PRETAIN=100,SVD=False):
    
    
    anom_eof = anom.copy()
    anom_eof[np.isnan(anom_eof)] = 0. #Put zeros in for the missing values.
    if SVD:

# Use SVD to get the EOFs. Transpose to do time x space.
        U,S,V = np.linalg.svd(anom_eof.T)
        pcs = U

        # Obtain EOFs by projecting the data matrix onto the PCs.               
        eofs = nancov(anom_eof,pcs)#.reshape((I,J,T),order='F') #EOFs with amplitude
        eigval = S**2/float(I*J) #Scale the singular values by N
    else:

# Use EIG to get the EOFs. Transpose to do time x space.
        Cov = nancov(anom_eof.T,anom_eof)/float(anom_eof.shape[0])
        eigval,B = np.linalg.eig(Cov)  # B = PCs, eigval = eigenvalues
        
        ###this should be in the else statement
        if eigval[-1] > eigval[1]: #Order the eigenvalues from largest to smallest.
            eigval = eigval[::-1]
            B = B[:,::-1]
        eofs = nancov(anom_eof,B)#.reshape((I,J,T),order='F')  #spatial EOF patterns
        pcs = B.copy()
# If PRETAIN is not 100, then choose the number of modes necessary
# to get to the percentage of the total variance you want to explain
# with the SVD analysis.
    eofs = eofs[:,:NMODES]
    pcs = pcs[:,:NMODES]
    eigval = eigval[:NMODES]
    vexp = 100.*eigval/np.sum(eigval)
    return eofs,pcs,vexp,eigval



#Call the compute EOF function to calculate the EOF (will want to change this for the concatonated matrix)
    
eof, pcs, vexp, eigval = Compute_EOFs(a.T, NMODES=16)


#Extract each field from multivariate EOF

zdreof = eof[0:39424]
kdpeof = eof[39424:78848]
precipeof = eof[78848:]


time,x,y = zdravg.shape
zdravg2 = zdravg.reshape((time,x*y), order = 'F').T
kdpavg2 = kdpavg.reshape((time,x*y), order = 'F').T

#%%

#Regress raw fields onto EOFs (regress precip_ano_no_nan)

stdPCs = pcs*np.nan; eofPatterns = precip_ano_no_nan*np.nan; stdECs = pcs*np.nan

I,J = lat.shape #precip latitude array
 
NMODES = 11
regressPatternPrecip = np.ones((I,J,NMODES))
tmpPrecip = np.ones((I*J,))*np.nan

for i in range(NMODES):
    stdPCs[:,i] = (pcs[:,i]-np.nanmean(pcs[:,i]))/np.nanstd(pcs[:,i])
    tmpPrecip[nanremove] = nancov(precip_ano_no_nan.T,stdPCs[:,i])
    regressPatternPrecip[:,:,i] = tmpPrecip.reshape(I,J)
    
    #eofPatterns[:,:,i] = regressPattern.reshape((I,J),order='F')
    I,J,T = regressPatternPrecip.shape
    precipdie = regressPatternPrecip.reshape((I*J,T), order = 'F')
    precipnandie = precipdie[nanremove]
    EC = nancov(precip_ano_no_nan,precipnandie)
    stdECs[:,i] = (EC[:,i]-np.nanmean(EC[:,i]))/np.nanstd(EC[:,i])

###Plot the EOF

def PlotEOF(pattern,lon,lat,mask,vexp,pcs,dates,REGION='Global',TIME='Monthly',MODE=1,DETREND=False):
    cmin = -1.; cmax = 1.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='BrBG',lut=nlevs)
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(211)
    print("LOOK:", pattern[:,:,MODE-1].shape)
    
#    m = Basemap(projection='npstere',boundinglat=35.,lon_0=-95.,llcrnrlat=lat[0,:].min(),\
#                urcrnrlat=lat[0,:].max(),llcrnrlon=lon[:,0].min(),\
#                urcrnrlon=lon[:,0].max(),resolution='l')
#    
    
    m = Basemap(projection='cyl',boundinglat=35.,lon_0=-95.,llcrnrlat=26,\
                urcrnrlat=33,llcrnrlon=-101,\
                urcrnrlon=-92,resolution='l')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x,y = m(lon,lat)
    print(x.shape, y.shape)
    cs = m.contourf(x,y,pattern[:,:,MODE-1],clevs,cmap=cmap,extend='both')
    #cmap._lut[nlevs/2-1:nlevs/2+1] = [1.,1.,1.,1.]  ###Centers colorbar at zero?
    cbar = m.colorbar(cs,pad='10%',location='bottom')
    cbar.ax.set_xlabel('Correlation',weight='bold',name='Calibri',size=14)
    cticks = []
    for i in clevs:
        cticks.append(int(i)) if i.is_integer() else cticks.append(i)
        cbar.set_ticks(clevs[::4])
        cbar.set_ticklabels(cticks[::4])
    for i in cbar.ax.yaxis.get_ticklabels():
        i.set_family('Calibri')
        i.set_size(14)
    ax.set_title('Precipitation (in) Mode 1 (15.7%)',name='Calibri',weight='bold',size=16)
    
    x2star,y2star = m(-97.3964,27.8006)
    m.plot(x2star,y2star,'ro',markersize=7, label = 'Corpus Christi')
    
    x3star,y3star = m(-95.3698,29.7604)
    m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
    plt.legend()
    
    #NH Winter SLP Anomalies Mode 1
    fig.tight_layout()
    plt.show(block=False) 

plotLon = lon
plotLat = lat
#plotMask = maskice[ilonice,:][:,ilatice]
plotSeries = pcs #if TIME == 'Monthly' else stdECs
plt.close('all')
for n in range(3):
    PlotEOF(regressPatternPrecip,plotLon,plotLat,vexp,pcs,precip_aug,plotSeries,MODE=n+1)





#ZDR

stdPCs = pcs*np.nan; eofPatterns = zdravg2*np.nan; stdECs = pcs*np.nan

I,J = radarlat2.shape #precip latitude array
 
NMODES = 11
regressPatternZdr = np.ones((I,J,NMODES))
tmpzdr = np.ones((I*J))*np.nan

for i in range(NMODES):
      
    stdPCs[:,i] = (pcs[:,i]-np.nanmean(pcs[:,i]))/np.nanstd(pcs[:,i])
    tmpzdr = nancov(zdravg2,stdPCs[:,i])
    regressPatternZdr[:,:,i] = tmpzdr.reshape(I,J)
    
    #eofPatterns[:,:,i] = regressPattern.reshape((I,J),order='F')
    I,J,T = regressPatternZdr.shape
    zdrdie = regressPatternZdr.reshape((I*J,T), order = 'F')
    #precipnandie = precipdie[nanremove]
    EC = nancov(zdravg2.T,zdrdie)
    stdECs[:,i] = (EC[:,i]-np.nanmean(EC[:,i]))/np.nanstd(EC[:,i])

###Plot the EOF

def PlotEOF(pattern,lon,lat,mask,vexp,pcs,dates,REGION='Global',TIME='Monthly',MODE=1,DETREND=False):
    cmin = -1.1; cmax = 1.1; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='bwr',lut=nlevs)
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(211)
    print("LOOK:", pattern[:,:,MODE-1].shape)
    
#    m = Basemap(projection='npstere',boundinglat=35.,lon_0=-95.,llcrnrlat=lat[0,:].min(),\
#                urcrnrlat=lat[0,:].max(),llcrnrlon=lon[:,0].min(),\
#                urcrnrlon=lon[:,0].max(),resolution='l')
#    
    
    m = Basemap(projection='cyl',boundinglat=35.,lon_0=-95.,llcrnrlat=26,\
                urcrnrlat=33,llcrnrlon=-101,\
                urcrnrlon=-92,resolution='l')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x,y = m(lon,lat)
    print(x.shape, y.shape)
    cs = m.contourf(x,y,pattern[:,:,MODE-1],clevs,cmap=cmap,extend='both')
    #cmap._lut[nlevs/2-1:nlevs/2+1] = [1.,1.,1.,1.]  ###Centers colorbar at zero?
    cbar = m.colorbar(cs,pad='10%',location='bottom')
    cbar.ax.set_xlabel('Correlation',weight='bold',name='Calibri',size=14)
    cticks = []
    for i in clevs:
        cticks.append(int(i)) if i.is_integer() else cticks.append(i)
        cbar.set_ticks(clevs[::4])
        cbar.set_ticklabels(cticks[::4])
    for i in cbar.ax.yaxis.get_ticklabels():
        i.set_family('Calibri')
        i.set_size(14)
    ax.set_title('3km $Z_{DR}$ (dB) Mode 1 (15.7%)',name='Calibri',weight='bold',size=16)
    x2star,y2star = m(-97.3964,27.8006)
    m.plot(x2star,y2star,'ro',markersize=7, label = 'Corpus Christi')
    
    x3star,y3star = m(-95.3698,29.7604)
    m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
    plt.legend()
    #NH Winter SLP Anomalies Mode 1
    fig.tight_layout()
    plt.show(block=False) 

plotLon = radarlon2
plotLat = radarlat2
#plotMask = maskice[ilonice,:][:,ilatice]
plotSeries = pcs #if TIME == 'Monthly' else stdECs
plt.close('all')
for n in range(3):
    PlotEOF(regressPatternZdr,plotLon,plotLat,vexp,pcs,precip_aug,plotSeries,MODE=n+1)




#KDP

stdPCs = pcs*np.nan; eofPatterns = kdpavg2*np.nan; stdECs = pcs*np.nan

I,J = radarlat2.shape #precip latitude array
 
NMODES = 10
regressPatternKdp = np.ones((I,J,NMODES))
tmpkdp = np.ones((I*J))*np.nan

for i in range(NMODES):
      
    stdPCs[:,i] = (pcs[:,i]-np.nanmean(pcs[:,i]))/np.nanstd(pcs[:,i])
    tmpkdp = nancov(kdpavg2,stdPCs[:,i])
    regressPatternKdp[:,:,i] = tmpkdp.reshape(I,J)
    
    #eofPatterns[:,:,i] = regressPattern.reshape((I,J),order='F')
    I,J,T = regressPatternKdp.shape
    kdpdie = regressPatternKdp.reshape((I*J,T), order = 'F')
    #precipnandie = precipdie[nanremove]
    EC = nancov(kdpavg2.T,kdpdie)
    stdECs[:,i] = (EC[:,i]-np.nanmean(EC[:,i]))/np.nanstd(EC[:,i])

###Plot the EOF

def PlotEOF(pattern,lon,lat,mask,vexp,pcs,dates,REGION='Global',TIME='Monthly',MODE=1,DETREND=False):
    cmin = -1.; cmax = 1.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
    nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='PRGn',lut=nlevs)
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(211)
    print("LOOK:", pattern[:,:,MODE-1].shape)
    
#    m = Basemap(projection='npstere',boundinglat=35.,lon_0=-95.,llcrnrlat=lat[0,:].min(),\
#                urcrnrlat=lat[0,:].max(),llcrnrlon=lon[:,0].min(),\
#                urcrnrlon=lon[:,0].max(),resolution='l')
#    
    
    m = Basemap(projection='cyl',boundinglat=35.,lon_0=-95.,llcrnrlat=26,\
                urcrnrlat=33,llcrnrlon=-101,\
                urcrnrlon=-92,resolution='l')
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    x,y = m(lon,lat)
    print(x.shape, y.shape)
    cs = m.contourf(x,y,pattern[:,:,MODE-1],clevs,cmap=cmap,extend='both')
    #cmap._lut[nlevs/2-1:nlevs/2+1] = [1.,1.,1.,1.]  ###Centers colorbar at zero?
    cbar = m.colorbar(cs,pad='10%',location='bottom')
    cbar.ax.set_xlabel('Correlation',weight='bold',name='Calibri',size=14)
    cticks = []
    for i in clevs:
        cticks.append(int(i)) if i.is_integer() else cticks.append(i)
        cbar.set_ticks(clevs[::4])
        cbar.set_ticklabels(cticks[::4])
    for i in cbar.ax.yaxis.get_ticklabels():
        i.set_family('Calibri')
        i.set_size(14)
    ax.set_title(' 3km $K_{DP}$ (deg/km) Mode 1 (15.7%)',name='Calibri',weight='bold',size=16)
    x2star,y2star = m(-97.3964,27.8006,)
    m.plot(x2star,y2star,'ro',markersize=7, label = 'Corpus Christi')
    
    x3star,y3star = m(-95.3698,29.7604)
    m.plot(x3star,y3star,'ko',markersize=7, label = 'Houston')
    plt.legend()
    #NH Winter SLP Anomalies Mode 1
    fig.tight_layout()
    plt.show(block=False) 

plotLon = radarlon2
plotLat = radarlat2
#plotMask = maskice[ilonice,:][:,ilatice]
plotSeries = pcs #if TIME == 'Monthly' else stdECs
plt.close('all')
for n in range(3):
    PlotEOF(regressPatternKdp,plotLon,plotLat,vexp,pcs,precip_aug,plotSeries,MODE=n+1)



#%%

def PlotEigValueSpectrum(eigval,N,NMODES=11,REGION='Global',TIME='Monthly'):

    #--------------------------------------------------------------------
    # Standard error in eigenvalue; 95% confidence interval; North et al. (1982)
    #--------------------------------------------------------------------
    deltaEIG = eigval[:NMODES]*np.sqrt(2./N)

    if eigval[:NMODES].max() > 50.:
        ylim = [0,40000.]; yticks = np.arange(0,400001.,5000.); #minorY = 5
    elif eigval[:NMODES].max() > 25:
        ylim = [0,50]; yticks = np.arange(0,51,10); #minorY = 2.
    else:
        ylim = [0,0.00003]; yticks = np.arange(0,0.00003,0.000005); #minorY = 2.

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)
    ax.set_position([0.1,0.2,0.5,0.5])
    ax.plot(np.arange(1,NMODES+1),eigval[:NMODES],'ro',lw=6)
    ax.errorbar(np.arange(1,NMODES+1),eigval[:NMODES],yerr=deltaEIG,color='r',fmt='o')
    ax.set(ylim=ylim,yticks=yticks,xlim=[0,NMODES+1],xticks=np.arange(1,NMODES+1))
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.tick_bottom()
    
    #ax.yaxis.set_minor_locator(MultipleLocator(minorY))
    ax.set_title('Eigenvalue Spectrum for the\n First Ten Modes', weight = 'bold', style = 'italic', size = 14)
    

    ax.set_ylabel(r'$\lambda$',name='Calibri',size=14,weight='bold')

    for i in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
        i.set_family('Calibri')
        i.set_size(14)

    plt.show(block=False)
 
    


N = 16
PlotEigValueSpectrum(eigval,N,NMODES=NMODES)

#%%

a = eigval[0:2]
b = eigval[2:]



#%%       
    

###Composite Difference between 3-hour precip, 3 hour averaged KDP, 3-hour averaged Zdr





      
        
        
        
        
        












    


    





    