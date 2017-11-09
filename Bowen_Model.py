# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 10:50:13 2017

@author: noahb
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 15:37:55 2017

@author: noahb
"""

import numpy as np
import matplotlib.pyplot as plt


r1 = 1*10**(-6) #initial drop size of radius in m
rho = 1000 #density of liquid water in units of kg/m^3
c = 1.0*10**(-3) # liquid water concentration in kg/m^3
c2 = 1.5*10**(-3)
w = 1 #vertical velocity in units of m/s
w2 = 2
e = 1 #collision efficiency
k = 8000 #in units per second

rstep = (0.002 - 0.000001)/100 #defining 100 values in our array for the final droplet size
r2 = np.arange(1*10**(-6), 2*10**(-3), rstep)  #final drop size array ranging from 1 micron to 2 mm


const1 = (4*rho)/(e*c)
const2 = w/k
const3 = (4*rho)/(e*c2)
const4 = w2/k
ratio = r2/r1

z = const1*(const2*np.log(ratio) - r2 + r1) #part b
z2 = const3*(const2*np.log(ratio) - r2 + r1) #part c
z3 = const1*(const4*np.log(ratio) - r2 + r1) #part d
z = z/1000
z2 = z2/1000
z3 = z3/1000



maxz = np.max(z) #find the maximum height
maxz2 = np.max(z2)
maxz3 = np.max(z3)


#maxheight = np.where(z==maxz)
#maxheight2 = np.where(z2==maxz2)
#maxheight3 = np.where(z3==maxz3)





###From the continuous-collection model
top = 4*rho*np.log(r2/r1)
bottom = k*c*1
bottom2 = k*c2*1
bottom3 = k*c*1

dt = (top/bottom)#/60
dt2 = (top/bottom2)#/60
dt3 = (top/bottom3)#/60

r2 = r2*10**(3) #r2 in units of mm

###change in time

####Plot######

FIG_WIDTH_INCHES = 10
FIG_HEIGHT_INCHES = 10

_, axes_object = plt.subplots(1, 1, figsize=(FIG_WIDTH_INCHES, FIG_HEIGHT_INCHES))
#plt.xscale('log')
xvals = dt
yvals = z
#plt.ylim(-1000, 2000)
plt.xlabel('Time (s)', size = 22.0)
plt.ylabel('Height (km)', size = 20.0)
axes_object.plot(xvals, yvals, 'b', linewidth = 5.0, label = 'Liquid water concentration = 1 g/^3')

xvals = dt2
yvals2 = z2
plt.xlabel('Time (s)', size = 22.0)
plt.ylabel('Height(km)', size = 20.0)
axes_object.plot(xvals, yvals2, 'r', linewidth = 5.0, label = 'Liquid water concentration = 1.5 g/m^3')


xvals = dt3
yvals3 = z3
plt.xlabel('Time (s)', size = 22.0)
plt.ylabel('Height (km)', size = 20.0)
axes_object.plot(xvals, yvals3, 'g', linewidth = 5.0, label = 'Vertical velocity = 2 m/s')


plt.legend()

plt.show()
