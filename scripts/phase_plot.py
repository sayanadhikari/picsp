import matplotlib.pyplot as plt
import numpy as np
import itertools
import seaborn as sns
from pylab import *

fig, ax = plt.subplots(2,2)


filename = "../output/phase_space/e0.dat"
i,xpos,xvel = np.loadtxt(filename, unpack=True)
x1 = xpos[i<30000]
v1 = xvel[i<30000]
x2 = xpos[i>30000]
v2 = xvel[i>30000]

xc = np.asarray([x1])
vc = np.asarray([v1])
xh = np.asarray([x2])
vh = np.asarray([v2])
ax[0,0].scatter(xc,vc,s=0.1,color='blue')
ax[0,0].scatter(xh,vh,s=0.1,color='red',alpha=0.5)
ax[0,0].set_title('timestep:0',fontsize = 12,loc='right')
ax[0,0].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
ax[0,0].legend(['cold','hot'])



filename = "../output/phase_space/e50.dat"
i,xpos,xvel = np.loadtxt(filename, unpack=True)

x1 = xpos[i<30000]
v1 = xvel[i<30000]
x2 = xpos[i>30000]
v2 = xvel[i>30000]

xc = np.asarray([x1])
vc = np.asarray([v1])
xh = np.asarray([x2])
vh = np.asarray([v2])
ax[0,1].scatter(xc,vc,s=0.1,color='blue')
ax[0,1].scatter(xh,vh,s=0.1,color='red',alpha=0.5)
ax[0,1].set_title('timestep:50',fontsize = 12,loc='right')
ax[0,1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))

filename = "../output/phase_space/e100.dat"
i,xpos,xvel = np.loadtxt(filename, unpack=True)
x1 = xpos[i<30000]
v1 = xvel[i<30000]
x2 = xpos[i>30000]
v2 = xvel[i>30000]

xc = np.asarray([x1])
vc = np.asarray([v1])
xh = np.asarray([x2])
vh = np.asarray([v2])
ax[1,0].scatter(xc,vc,s=0.1,color='blue')
ax[1,0].scatter(xh,vh,s=0.1,color='red',alpha=0.5)
ax[1,0].set_title('timestep:100',fontsize = 12,loc='right')
ax[1,0].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
ax[1,0].legend(['cold','hot'])



filename = "../output/phase_space/e150.dat"
i,xpos,xvel = np.loadtxt(filename, unpack=True)

x1 = xpos[i<30000]
v1 = xvel[i<30000]
x2 = xpos[i>30000]
v2 = xvel[i>30000]

xc = np.asarray([x1])
vc = np.asarray([v1])
xh = np.asarray([x2])
vh = np.asarray([v2])
ax[1,1].scatter(xc,vc,s=0.1,color='blue')
ax[1,1].scatter(xh,vh,s=0.1,color='red',alpha=0.5)
ax[1,1].set_title('timestep:150',fontsize = 12,loc='right')
ax[1,1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))



plt.show()
