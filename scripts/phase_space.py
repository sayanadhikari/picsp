#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin

fig, ax = plt.subplots(3,1)

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
xe0 = []
vxe0 = []
ke0 = []
xpos = []
xvel = []

for i in range(0,1000,5):
    datae = h5["particle.e/%d"%i]
    xe = datae[:,0]
    ye = datae[:,1]
    vxe = datae[:,2]
    vye = datae[:,3]
    x0 = datae[0,0]
    vx0 = datae[0,2]
    # xpos.append(xe)
    # xvel.append(vxe)
    xe0.append(x0)
    vxe0.append(vx0)
    ke0.append(vx0**2)
    KE.append(np.sum(vxe**2) + np.sum(vye**2))

KE = np.array(KE)
time = np.arange(0,1000,5)
position = np.arange(0,0.784,33)
xe0 = np.array(xe0)
vxe0 = np.array(vxe0)
ke0 = np.array(ke0)
# xpos = np.array(xpos)
# xvel = np.array(xvel)
# print(KE)
# print(xvel)
# plt.scatter(xpos,xvel,s=0.5)

ax[0].plot(time,xe0)
ax[0].set_title("x vs t")
ax[1].plot(time,ke0)
ax[1].set_title("ke vs t")
ax[2].plot(xe0,vxe0)

# plt.plot(xe0,vxe0)

# Phi = []
# Efx = []
# # POT = []
# j=50
# # for j in range(0,2000,5):
# data_phi = h5["phi/%d"%j]
# data_efx = h5["efx/%d"%j]
# phi = data_phi[0,1]
# efx = data_efx[:]
# Phi.append(phi)
# Efx.append(efx)
# Efx = np.array(Efx)
# POT.append(np.sum(phi*rho))
#
# Phi = np.array(Phi)
# Rho = np.array(rho)
# POT = np.array(POT)
#ax[1,1].plot(time,Phi)
# print(Efx)
# plt.plot(position,efx)
# plt.legend()
# plt.show()

#plt.plot(xe0,vxe0)
plt.show()
