# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 10:43:25 2017

@author: noahb
"""

import numpy as np
import matplotlib.pyplot as plt




g = 9.8
h = 10**4 #in meters
c = np.sqrt(g*h)
vel = 7.292*10**(-5) ##per second
re = 6371000 # in meters
beta = (2*vel)/re
wavenum = np.arange(0, (np.pi*(1*10**-6))/2, 1*10**-7 )


mu = wavenum*((c/beta)**(1/2))



def calcnu(omega, beta, c):
    nu = omega/(np.sqrt(beta*c))
    return(nu)
    

####Poincare waves####
def poincare(c, wavenum, beta, n):
    freq1 = np.sqrt((c**2)*(wavenum**2) + (2*n +1)*beta*c )
    return(freq1)

def poincare2(c, wavenum, beta, n):
    freq1 = -np.sqrt((c**2)*(wavenum**2) + (2*n +1)*beta*c )
    return(freq1)
        
    

#### Kelvin Waves ####
def kelvin(wavenum, c):
    freq2 = wavenum*c
    return(freq2)

#### Mixed ####
##Positive 
def mixedplus(c, wavenum, beta):
    freq3 = ((wavenum*c)/2) +((((wavenum**2)*(c**2))) + beta*c)**(1/2)
    return(freq3)
##Negative
def mixedneg(c, wavenum, beta):
    freq4 = ((wavenum*c)/2) + -1 *((((wavenum**2)*(c**2))) + beta*c)**(1/2)
    return(freq4)

#### Equatorial Rossby Waves ####
def eqros(c, wavenum, beta, n): 
    freq5 = (-1*beta*wavenum)/(((wavenum**2)+(2*n + 1)*(beta/c)))
    return(freq5)


    
#####plot############

FIG_WIDTH_INCHES = 10
FIG_HEIGHT_INCHES = 10

_, axes_object = plt.subplots(1, 1, figsize=(FIG_WIDTH_INCHES, FIG_HEIGHT_INCHES))

plt.xlim(0,4)
plt.ylim(-5,5)
for n in[1, 2, 3]:
    freq1 = poincare(c, wavenum, beta, n)
    nu = calcnu(freq1, beta, c)
    plt.plot(mu, nu, 'r-', label = 'Poincare Wave')
    freq1 = poincare2(c, wavenum, beta, n)
    nu = calcnu(freq1, beta, c)
    plt.plot(mu, nu, 'r-', label = 'Poincare Wave')
    
freq2 = kelvin(wavenum, c)
nu = calcnu(freq2, beta, c)
plt.plot(mu, nu, 'b-', label = 'Kelvin Wave')
freq3 = mixedplus(c, wavenum, beta)
nu = calcnu(freq3, beta, c)
plt.plot(mu, nu, 'g-', label = 'Mixed Rossby Gravity Waves')
freq4 = mixedneg(c, wavenum, beta)
nu = calcnu(freq4, beta, c)
plt.plot(mu, nu, 'g-', label = 'Mixed Rossby Gravity Waves')

for n in[1, 2, 3]:
    freq5 = eqros(c, wavenum, beta, n)
    nu = calcnu(freq5, beta, c)
    plt.plot(mu, nu, 'p-', label = 'Equatorial Rossby Waves')
plt.xlabel('mu', size = 16)
plt.ylabel('nu', size = 16)
plt.legend()
plt.show()

