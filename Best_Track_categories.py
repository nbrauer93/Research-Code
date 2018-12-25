# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 23:18:16 2018

@author: noahb
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap,maskoceans,interp
from netCDF4 import Dataset, num2date, MFDataset




track = 'Harvey_Best_Track.txt'

harvey = np.loadtxt(track, dtype = 'str', skiprows = 1)

#Lat-lon of actual hurricane location

latitude = harvey[:,4] 
longitude = harvey[:,5]
time = harvey[:,:2]
pres = harvey[:,7]
wind = harvey[:,6]

#Meshgrid lat and lon 

lat2,lon2 = np.meshgrid(latitude,longitude)[:]

wind_float = np.ones(74)*np.nan

for i in range(len(wind)):
    string = wind[i]
    strLen = len(string)
    wind_float[i] = string[0:strLen - 1]
    
    

#Convert from degrees to km

def latlon_to_dist(lat_input2,lat_input1,lon_input2,lon_input1):
    distlat = (lat_input2-lat_input1)*111
    distlon = (lon_input2 - lon_input1)*111
    totaldist = np.sqrt((distlat)**2 + (distlon)**2)
    
    return(totaldist)
        

#These are the coordinates of Harvey
    
coords = np.array([latitude,longitude],dtype = 'float').T # 74x2

#Designate categories for TS, CAT 1, CAT 2, etc.
trop_dep = np.where(wind_float<39)[0]
trop_storm = np.where((wind_float>38)&(wind_float<=74))[0]
category_1 = np.where((wind_float>74)&(wind_float<=95))[0]
category_2 = np.where((wind_float>96)&(wind_float<=110))[0]
category_3 = np.where((wind_float>110)&(wind_float<=129))[0]
category_4 = np.where((wind_float>129)&(wind_float<=156))[0]

#Assign these to the wind array

td = wind_float[trop_dep]
ts = wind_float[trop_storm]
cat1 = wind_float[category_1]
cat2 = wind_float[category_2]
cat3 = wind_float[category_3]
cat4 = wind_float[category_4]

#Assign correct lat-lons to each category

td_lat,td_lon = latitude[trop_dep], longitude[trop_dep]
ts_lat,ts_lon = latitude[trop_storm], longitude[trop_storm]
cat1_lat,cat1_lon = latitude[category_1], longitude[category_1]
cat2_lat,cat2_lon = latitude[category_2], longitude[category_2]
cat3_lat,cat3_lon = latitude[category_3], longitude[category_3]
cat4_lat,cat4_lon = latitude[category_4], longitude[category_4]


#Define arrays for each

td_coords = np.array([td_lat,td_lon], dtype = 'float').T
ts_coords = np.array([ts_lat,ts_lon], dtype = 'float').T
cat1_coords = np.array([cat1_lat,cat1_lon], dtype = 'float').T
cat2_coords = np.array([cat2_lat,cat2_lon], dtype = 'float').T
cat3_coords = np.array([cat3_lat,cat3_lon], dtype = 'float').T
cat4_coords = np.array([cat4_lat,cat4_lon], dtype = 'float').T

#%%

plt.figure(figsize=(10,10))

#Largeest domain (entire lifetime)
#xlim = np.array([-105,-40]); ylim = np.array([12,40])

#Uncomment for smaller domain
#xlim = np.array([-105,-40]); ylim = np.array([15,40])

#Local domaian
xlim = np.array([-105,-85]); ylim = np.array([20,40])


m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='l')
m.drawcoastlines(); m.drawstates(), m.drawcountries()  
cs = m(lon2,lat2) 
plt.title('Hurricane Harvey Track', size = 22, weight = 'bold')


tdjunk_x,tdjunk_y = m(td_coords[:,1],td_coords[:,0])
m.plot(tdjunk_x,tdjunk_y,'ko',markersize = 8, label = 'TD')

tsjunk_x,tsjunk_y = m(ts_coords[:,1],ts_coords[:,0])
m.plot(tsjunk_x,tsjunk_y,'bo',markersize = 8, label = 'TS')

cat1junk_x,cat1junk_y = m(cat1_coords[:,1],cat1_coords[:,0])
m.plot(cat1junk_x,cat1junk_y,'go',markersize = 8, label = 'Cat 1')

cat2junk_x,cat2junk_y = m(cat2_coords[:,1],cat2_coords[:,0])
m.plot(cat2junk_x,cat2junk_y,'yo',markersize = 8, label = 'Cat 2')


cat3junk_x,cat3junk_y = m(cat3_coords[:,1],cat3_coords[:,0])
m.plot(cat1junk_x,cat1junk_y,'mo',markersize = 8, label = 'Cat 3')

'''
cat4junk_x,cat4junk_y = m(cat4_coords[:,1],cat4_coords[:,0])
m.plot(cat4junk_x,cat4junk_y,'ro',markersize = 8, label = 'Cat 4')
'''


plt.legend()
plt.show()






