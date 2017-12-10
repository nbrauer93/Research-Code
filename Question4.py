
import matplotlib.pyplot as plt
import numpy as np


#define the courant number



mu = np.linspace(-5, 5, 100)
#define a spatial resolution from 0 to pi
kdeltax = np.linspace(0, np.pi, 100)

#ampfactor = 1 - 2*mu*(1 - mu)*(1 - np.cos(kdeltax))

dummy = np.zeros((100, 100))
for i in range (0,100):
    for j in range (0,100):
        dummy[i, j] = np.sqrt(((mu[i]**2)*((np.sin(kdeltax[j]))**2))+((1-((mu[i]**2)*(1-(np.cos(kdeltax[j])))))**2))
        
        
ampfactor = dummy
       



#####plot############

FIG_WIDTH_INCHES = 10
FIG_HEIGHT_INCHES = 10

_, axes_object = plt.subplots(1, 1, figsize=(FIG_WIDTH_INCHES, FIG_HEIGHT_INCHES))

plt.ylabel('Courant Number')
plt.xlabel('KDeltax')

y1 = ampfactor
x , y = np.meshgrid(kdeltax, mu)
level = np.arange(0, 2, 0.01) 
plt.contourf(x, y, y1, level, cmap = 'coolwarm', extend = 'both')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Amplification Factor')
legend = plt.legend()
plt.title('Lax-Wendroff')
plt.show()

