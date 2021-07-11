#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin

# fig, axs = plt.subplots(3,1)
# axs = axs.ravel()

path = pjoin("..","output")
print(path)

file_name = 'data'

h5 = h5py.File(pjoin(path,file_name+'.h5'),'r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
# Lz = h5.attrs["Lz"]

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

KE = []
xpos = []
xvel = []

# for j in range(1,3):
for i in range(0,1000,500):
    for j in range(1,3):
        fig, axs = plt.subplots(j,1)
        datae = h5["particle.e/%d"%i]
        xe = datae[:,0]
        ye = datae[:,1]
        vxe = datae[:,2]
        vye = datae[:,3]
        xpos.append(xe)
        xvel.append(vxe)
        KE.append(np.sum(vxe**2) + np.sum(vye**2))
        axs[j-1].scatter(xpos,xvel,s=0.5)


# i = 500
# datae = h5["particle.e/%d"%i]
# xe = datae[:,0]
# vxe = datae[:,2]
# xpos.append(xe)
# xvel.append(vxe)
# plt.scatter(xpos,xvel,s=0.05)


# i = 400
# datae = h5["particle.e/%d"%i]
# xe = datae[:,0]
# ye = datae[:,1]
# vxe = datae[:,2]
# vye = datae[:,3]
# xpos.append(xe)
# xvel.append(vxe)
# ax[0,1].scatter(xpos,xvel,s=0.1)
#
#
# i = 500
# datae = h5["particle.e/%d"%i]
# xe = datae[:,0]
# ye = datae[:,1]
# vxe = datae[:,2]
# vye = datae[:,3]
# xpos.append(xe)
# xvel.append(vxe)
# ax[1,0].scatter(xpos,xvel,s=0.1)
#
#
# i = 600
# datae = h5["particle.e/%d"%i]
# xe = datae[:,0]
# ye = datae[:,1]
# vxe = datae[:,2]
# vye = datae[:,3]
# xpos.append(xe)
# xvel.append(vxe)
# ax[1,1].scatter(xpos,xvel,s=0.1)


# plt.legend()
plt.show()

#plt.plot(xe0,vxe0)
#plt.show()
