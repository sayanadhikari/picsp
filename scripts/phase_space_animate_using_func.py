import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from pylab import *

nTimeSteps = 1000


DATA = []
# DEFAULT VALUES
eps_0   = 8.85418782e-12
q       = 1.602e-19 # charge
M_i     = 1.673e-27; #mass Ions
M_e     = 9.109e-31#4e-29
Ti      = 0.026
Te      = 1
TeSI    = sqrt(2*(Te*q)/M_e)
TiSI    = sqrt(2*(Ti*q)/M_i)
gamma   = 1.67
KB      = 1.38064852E-23
n_0     = 1e16

DL       = np.sqrt((TeSI*KB*eps_0)/(n_0*q*q))

# numerical data file
DIR ="../output/phase_space"

for i in range(0,nTimeSteps,50):
    filename1=DIR+'/i%d'%i+'.dat'
    filename2=DIR+'/e%d'%i+'.dat'

    x1,v1 = np.loadtxt(filename1, unpack=True)
    x2,v2 = np.loadtxt(filename2, unpack=True)
 
 DATA = DATA[start_index:,:,:]
 ani=animation.FuncAnimation(fig,animate,len(DATA[:,0,0]),interval=interval*1e+3,blit=False)
 
 
figure(1)
plt.scatter(x1/DL,v1,s=1,marker='.')
plt.xlabel("$x/\\lambda_D$")
plt.ylabel("$V_x$")
plt.title('Ion Phase Space')
plt.tight_layout()
plt.pause(0.1)
plt.clf()


figure(2)
plt.scatter(x2/DL,v2,s=1,marker='.')
plt.xlabel("$x/\\lambda_D$")
plt.ylabel("$V_x$")
plt.title('Electron Phase Space')
plt.tight_layout()
plt.pause(0.1)
plt.clf()

