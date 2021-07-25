import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pylab import *
import h5py
from os.path import join as pjoin

fig, ax = plt.subplots()


path = pjoin("..","output")
print(path)

file_name = 'data'

h5 = h5py.File(pjoin(path,file_name+'.h5'),'r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
# Lz = h5.attrs["Lz"]

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

xe0 = []
vxe0 = []

i=70
# for i in range(0,200,1):
datae = h5["particle.e/%d"%i]
xe = datae[:,0]
ye = datae[:,1]
vxe = datae[:,2]
vye = datae[:,3]
xe0.append(xe)
vxe0.append(vxe)

ax.plot(xe0,vxe0,'.')
plt.show()

# ax[0,0].plot(xpos,xvel,'o',markersize=0.1)
# ax[0,0].set_ylabel('v/$v_0$',fontsize = 8)
# ax[0,0].set_title('t_160',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e50.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[0,1].plot(xpos,xvel,'o',markersize=0.1)
# ax[0,1].set_title('t_50',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e100.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[0,2].plot(xpos,xvel,'o',markersize=0.1)
# ax[0,2].set_title('t_100',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e150.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[0,3].plot(xpos,xvel,'o',markersize=0.1)
# ax[0,3].set_title('t_150',fontsize = 8,loc='right')
#
# filename = "output/phase_space/e200.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[1,0].plot(xpos,xvel,'o',markersize=0.1)
# ax[1,0].set_ylabel('v/$v_0$',fontsize = 8)
# ax[1,0].set_title('t_200',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e250.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[1,1].plot(xpos,xvel,'o',markersize=0.1)
# ax[1,1].set_title('t_250',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e500.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[1,2].plot(xpos,xvel,'o',markersize=0.1)
# ax[1,2].set_title('t_500',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e1000.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[1,3].plot(xpos,xvel,'o',markersize=0.1)
# ax[1,3].set_title('t_1000',fontsize = 8,loc='right')
# ax[1,3].set_ylim([10,-10])
#
# filename = "output/phase_space/e1500.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[2,0].plot(xpos,xvel,'o',markersize=0.1)
# ax[2,0].set_ylabel('v/$v_0$',fontsize = 8)
# ax[2,0].set_title('t_1500',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e2000.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[2,1].plot(xpos,xvel,'o',markersize=0.1)
# ax[2,1].set_title('t_2000',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e2500.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[2,2].plot(xpos,xvel,'o',markersize=0.1)
# ax[2,2].set_title('t_2500',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e3000.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[2,3].plot(xpos,xvel,'o',markersize=0.1)
# ax[2,3].set_title('t_3000',fontsize = 8,loc='right')
#
#
# filename = "output/phase_space/e3500.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[3,0].plot(xpos,xvel,'o',markersize=0.1)
# ax[3,0].set_ylabel('v/$v_0$',fontsize = 8)
# ax[3,0].set_xlabel('2$\pi$x/L',fontsize = 8)
# ax[3,0].set_title('t_3500',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e4000.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[3,1].plot(xpos,xvel,'o',markersize=0.1)
# ax[3,1].set_xlabel('2$\pi$x/L',fontsize = 8)
# ax[3,1].set_title('t_4000',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e4500.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[3,2].plot(xpos,xvel,'o',markersize=0.1)
# ax[3,2].set_xlabel('2$\pi$x/L',fontsize = 8)
# ax[3,2].set_title('t_4500',fontsize = 8,loc='right')
#
#
#
# filename = "output/phase_space/e5000.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename, unpack=True)
# ax[3,3].plot(xpos,xvel,'o',markersize=0.1)
# ax[3,3].set_xlabel('2$\pi$x/L',fontsize = 8)
# ax[3,3].set_title('t_5000',fontsize = 8,loc='right')
# ax[3,3].set_ylim([10,-10])


# """
# filename2 = "output/phase_space/e10.dat"
# xpos,ypos,xvel,yvel = np.loadtxt(filename2, unpack=True)
# subplot(442)
# plt.plot(xpos,xvel,'o',markersize=1),
# plt.title('t_25',fontsize = 8,loc='right')
# plt.xticks(fontsize=5)
# plt.yticks(fontsize=5)
# plt.ylabel('v/$v_0$',fontsize = 8)
# plt.xlabel('2$\pi$x/L',fontsize = 8)
# #plt.ylim([-2, 2])
# """

plt.show()
