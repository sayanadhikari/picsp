#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pylab import *
# numerical data file
DIR ="../output/"

file_name = "data"#"rhoNeutral" #"P"

h5 = h5py.File(DIR+file_name+'.h5','r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]

Nx = int(h5.attrs["Nx"])
Ny = int(h5.attrs["Ny"])

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
# X, Y = np.meshgrid(x, y)

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

data =  h5["phi/%d"%Nt]
datae = h5["particle.e/%d"%data_num[i]]

xLen,iDen,eDen,iVel,eVel,rho,phi,ef = np.loadtxt(filename, unpack=True)

ActiveSubplot = True


if nTimePhase == 0:
    inix = nTimePhase*numCells
inix = nTimePhase*numCells+nTimePhase


if ActiveSubplot==False:
    figure(1)
    plt.plot(xLen[inix:inix+numCells],iDen[inix:inix+numCells], linewidth=1,color='blue')
    plt.plot(xLen[inix:inix+numCells],eDen[inix:inix+numCells], linewidth=1,color='red')
    plt.legend(('$N_i$', '$N_e$'))
    # plt.plot(xLen,eDen, linewidth=4,color='blue')
    plt.xlabel('Length (m)')
    plt.ylabel('Density ($m^{-3}$)')
    plt.title('Electron / Ion Density')

    figure(2)
    plt.plot(xLen[inix:inix+numCells],ef[inix:inix+numCells], linewidth=1,color='red')
    plt.xlabel('Length (m)')
    plt.ylabel('$E$ (V/m)')
    plt.title('Electric Field')

    figure(3)
    plt.plot(xLen[inix:inix+numCells],rho[inix:inix+numCells], linewidth=1,color='red')
    plt.xlabel('Length (m)')
    plt.ylabel('$\\rho$ (V/m)')
    plt.title('Charge Density')

    figure(4)
    plt.plot(xLen[inix:inix+numCells],phi[inix:inix+numCells], linewidth=1,color='blue')
    plt.xlabel('Length (m)')
    plt.ylabel('$\\phi$ (V)')
    plt.title('Potential')


if ActiveSubplot==True:
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.subplot(221)
    plt.plot(xLen[inix:inix+numCells],iDen[inix:inix+numCells], linewidth=1,color='blue')
    plt.plot(xLen[inix:inix+numCells],eDen[inix:inix+numCells], linewidth=1,color='red')
    plt.legend(('$N_i$', '$N_e$'))
    # plt.plot(xLen,eDen, linewidth=4,color='blue')
    plt.xlabel('Length (m)')
    plt.ylabel('Density ($m^{-3}$)')
    plt.title('Electron / Ion Density')
    plt.subplot(222)
    plt.plot(xLen[inix:inix+numCells],phi[inix:inix+numCells], linewidth=1,color='blue')
    plt.xlabel('Length (m)')
    plt.ylabel('$\\phi$ (V)')
    plt.title('Potential')
    plt.subplot(223)
    plt.plot(xLen[inix:inix+numCells],ef[inix:inix+numCells], linewidth=1,color='red')
    plt.xlabel('Length (m)')
    plt.ylabel('$E$ (V/m)')
    plt.title('Electric Field')
    plt.subplot(224)
    plt.plot(xLen[inix:inix+numCells],rho[inix:inix+numCells], linewidth=1,color='green')
    plt.xlabel('Length (m)')
    plt.ylabel('$\\rho$ (V/m)')
    plt.title('Charge Density')

plt.show()
