# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 18:37:29 2017

@author: noahb
"""

import numpy as np

#Define range of temperature

tupper = 273.15 - 5
tlower = 273.15 - 25

#Define equal increments of temperature

tempincrement = (tupper - tlower)/1000 

#Create an array of temperatures within given range by the assigned increment (tempincrement)
temp = np.arange(tlower, tupper, tempincrement)


#Do the same for z as we did with temperature
zupper = 10
zlower = 0
#zincrement = (zlower - zupper)/1000   ##python doesn't like this so I manually entered in 10/1000 into the next line
z = np.arange(zlower, zupper, 0.01)
#Created an array of all z values with 1000 increments ranging from 0 to 10 cm



#Do the messy arithmetic in steps
t1w = -7.90298*((373.15/temp) - 1)
#print (t1w)
t2w = 5.02808
t3w = np.log10(373.15/temp)
#print(t3w)
t4w = -1.3816*10**(-7)
t5w = temp/373.15
t6w = 1 - t5w
t7w = (10**(11.344*t6w) - 1)
#print(t7w)
t8w = 8.1328*10**-3
t9w = (373.15/temp) - 1
t10w = (10**(-3.49148*t9w))- 1
#print(t10w)
t11w = np.log10(1013.25)


## Saturation vapor pressure for water 

def essat(temp):
    sat_vap_water = t1w + t2w*t3w - t4w*t7w + t8w*t10w + t11w
    sat_vap_water = 10**(sat_vap_water)
    #print("Saturation vapor pressure for water is", sat_vap_water, "hPa")
    
    return sat_vap_water

watersat = essat(temp)
#return the saturation vapor pressure for water --> Use watersat for this variable



## Saturation vapor pressure for ice

#Messy Groff-Gratch arithmetic 
t1i = -9.09718*((273.15/temp) - 1)
t2i = -3.56654*np.log10(273.15/temp)
t3i = 0.876793*(1 - (temp/273.15))
t4i = np.log10(6.11)

def eisat(temp):
    sat_vap_ice = t1i + t2i + t3i + t4i
    sat_vap_ice = 10**(sat_vap_ice)
    #print("Saturation vapor pressure for ice is", sat_vap_ice, "hPa")
    return sat_vap_ice

icesat = eisat(temp)

#return the saturation vapor pressure for ice --> Use icesat for this variable


##### Compute the maximum supersaturation in the chamber
#Assume vapor pressure (e) varies linearly with height

#define each step (interval) for icesat, with 1000 steps for vapor pressure of ice
intervalice = (icesat[-1]-icesat[0])/1000

#define each step (interval) for watersat, with 1000 steps for vapor pressure of water
intervalwater = (watersat[-1]- watersat[0])/1000



#create a linear function for the variation of vapor pressure of ice with respect to height
vappice = np.arange(icesat[0], icesat[-1], intervalice)

#create a linear function for the variation of vapor pressure of water with respect to height
vappwater = np.arange(watersat[0], watersat[-1], intervalwater)

maxSSwater = (vappwater/watersat) - 1 #Supersaturation equation for water
maxSSice = (vappice/icesat) - 1 #Supersaturation equation for ice

#Calculate the maximum saturation vapor pressure values for water and ice

maxsice = np.max(maxSSice)  
maxswater = np.max(maxSSwater)

print(maxsice, "hPa") #Max supersaturation of ice
print(maxswater, "hPa") #Max supersaturation of water


max_height_ice = np.where(maxSSice==maxsice)
zmaxice = z[max_height_ice]  #applying the maximum supersaturation value of ice to the z(height) array that we created

max_height_water = np.where(maxSSwater==maxswater)
zmaxwater = z[max_height_water]  #applying the maximum supersaturation value of water to the z(height) array that we created

print("The height of maximum supersaturation for ice is", zmaxice, "cm")
print("The height of maximum supersaturation for water is", zmaxwater, "cm")










    
