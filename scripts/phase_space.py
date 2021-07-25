#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt
from os.path import join as pjoin

fig, ax1 = plt.subplots(2,2)

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

for i in range(0,2500,5):
    datae = h5["particle.e/%d"%i]
    xe = datae[:,0]
    ye = datae[:,1]
    vxe = datae[:,2]
    vye = datae[:,3]
    x0 = datae[4,0]
    vx0 = datae[4,2]
    # xpos.append(xe)
    # xvel.append(vxe)
    xe0.append(x0)
    vxe0.append(vx0)
    ke0.append(vx0**2)
    KE.append(np.sum(vxe**2) + np.sum(vye**2))

KE = np.array(KE)
time = np.arange(0,2500,5)
position = np.linspace(0,1,65)
position = np.reshape(position,(65,1))
xe0 = np.array(xe0)
vxe0 = np.array(vxe0)
ke0 = np.array(ke0)
# xpos = np.array(xpos)
# xvel = np.array(xvel)
# print(KE)
# print(xvel)
# plt.scatter(xpos,xvel,s=0.5)
# plot1 = plt.figure(1)
# plt.plot(time,ke0)

# ax1[0,0].plot(time,xe0)
# ax1[0,0].set_title("x vs t (Particle No. 4)")
# ax1[0,1].plot(time,ke0)
# ax1[0,1].set_title("KE (Particle No. 4)")
# ax1[1,0].plot(time,KE)
# ax1[1,0].set_title("Total KE")
# ax1[1,1].plot(xe0,vxe0)

# plt.plot(xe0,vxe0)
Rho = []
Phi = []
Efx = []
j=2500
# for j in range(0,2500,5):
data_phi = h5["phi/%d"%j]
data_efx = h5["efx/%d"%j]
data_rho = h5["rho/%d"%j]
phi = data_phi[:,1]
efx = data_efx[:,1]
rho = data_rho[:,1]
Phi.append(phi)
Efx.append(efx)
Rho.append(rho)
Phi = np.array(Phi)
Phi = np.ravel(Phi)
Efx = np.array(Efx)
Efx = np.ravel(Efx)
Rho = np.array(Rho)
Rho = np.ravel(Rho)
# Phi = np.reshape(Phi,(17,3))
# pot = np.array(Phi)*np.array(Rho)
# plot2 = plt.figure(2)

ax1[0,0].plot(position, Phi)
ax1[0,0].set_title("x vs Phi")
ax1[0,1].plot(position,Rho)
ax1[0,1].set_title("x vs Rho")
ax1[1,0].plot(position, Efx)
ax1[1,0].set_title("x vs Efx")
ax1[1,1].plot(xe0,vxe0)


# plt.plot(position,Efx)
# plt.xlabel("x")
# plt.ylabel("Efx")
# plt.title("Electric field at ts = 100")
# print(Phi)
plt.grid()
plt.show()
