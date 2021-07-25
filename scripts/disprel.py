#!/usr/bin/env python3

# Usage: ./disprel.py path/to/folder

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mp
import numpy as np
import sys
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
from tqdm import tqdm
import argparse
import h5py




####################################################





folder = pjoin('..','output')

parser = argparse.ArgumentParser(description='Plasma Dispersion Processor')
parser.add_argument('-per','--periodic', action='store_true', help='Add this if the system is periodic in Y')
parser.add_argument('-yLoc','--yLocation', default=1, type=int, help='In bounded (in Y) system Choose Y location, Options: e.g. 1 (Any number between 0-Ny)')
parser.add_argument('-pl','--plot', action='store_true', help='Add this if you want to plot the figure')
parser.add_argument('-n','--norm', default='omega_pe', type=str, help='Normalizing frequency, Options: omega_pi, omega_pe')


args        = parser.parse_args()
periodic    = args.periodic
yLoc        = args.yLocation	#8 # Choose y location
plot        = args.plot
norm        = args.norm

# Set processed data directory
# folder_base= os.path.basename(os.path.dirname(folder))
savedir     = pjoin(folder, 'processed_Efield')

h5 = h5py.File(pjoin(folder,'data.h5'),'r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]

Nx = int(h5.attrs["Nx"])
Ny = int(h5.attrs["Ny"])


dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

ne = h5.attrs["density"]
ni = ne

mi  = h5.attrs["massI"]

Te        = h5.attrs["vthE"]             # in eV
Ti        = h5.attrs["vthI"]           # in eV
units = 'EV'
#######################################################
# Analytical ion-acoustic dispresion relation


eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
e = constants('elementary charge')
c0 = constants('speed of light in vacuum')


gamma_e = 5./3


if units == 'MKS':
  tEeV  = 0.5*me*(Te*Te)/e
  tEK   = tEeV*11604.525
  tIeV  = 0.5*mi*(Ti*Ti)/e
  tIK   = tIeV*11604.525
if units == 'EV':
  tEK    = Te*11604.525
  tESI   = np.sqrt(2*(Te*e)/me)
  tIK   = Ti*11604.525
  tISI   = np.sqrt(2*(Te*e)/mi)

Te  = tEK #vars['tEK'] #1.6*11604
Ti  = tIK #vars['tIK'] #0.1*11604

wpi = np.sqrt(e**2*ni/(eps0*mi))
wpe = np.sqrt(e**2*ne/(eps0*me))
dl	= np.sqrt(eps0*kb*Te/(ni*e*e))
dli	= np.sqrt(eps0*kb*Ti/(ni*e*e))
print("wpe={}".format(wpe))
print("dl={}".format(dl))
cia = np.sqrt(gamma_e*kb*Te/mi)





x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X, Y = np.meshgrid(x, y)

dt = h5.attrs["dt"]*dp

n = np.arange(start=0, stop=Nt, step=dp, dtype=int)
t = n*dt


if  os.path.exists(savedir) and os.path.exists(pjoin(savedir,'pro_data_file.npz')):
        print('processed data exists. Loading data ...')
        data = np.load(pjoin(savedir,'pro_data_file.npz'))
        f = data['data']
        x = data['x']
        print('Shape of loaded data fp: ',f.shape)
        print("Shape of loaded data x",x.shape)

else:
  if  os.path.exists(savedir)== False:
    os.mkdir(savedir)
  print('processed data does not exist. Wait, processing data ...')

  # dataset index
  data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

  # READ FIELD
  f = np.zeros((Nt, Nx, Ny))

  for i in tqdm(range(len(data_num))):
      f[i,:,:] = h5['efx/%d'%data_num[i]]
  # exit()
    # Remove y
  if periodic == False:
      f = f[:,:,yLoc-10:yLoc+10]
  f = np.average(f, axis=2)
  # x = np.average(x, axis=1)
  np.savez_compressed(pjoin(savedir,'pro_data_file.npz'),data=f,x=x)

if plot:
    # FFT axes
  # dt = vars['timeStep']*vars['save_step'] #t[1]-t[0]


  dx = Lx/Nx #x[1]-x[0]
  omega = 2*np.pi*np.arange(Nt)/(Nt*dt)
  k     = 2*np.pi*np.arange(Nx)/(Nx*dx)
  print('Length of k: ',len(k))
  print('Max of k: ',np.max(k))
  Omega, K = np.meshgrid(omega, k, indexing='ij')
  print('Shape of Omega: ',Omega.shape)
  F = np.fft.fftn(f, norm='ortho')

  halflen = np.array(F.shape, dtype=int)//2
  Omega = Omega[:halflen[0],:halflen[1]]
  K = K[:halflen[0],:halflen[1]]
  F = F[:halflen[0],:halflen[1]]
  nK  = Nx

  ka = np.linspace(0, np.max(K), nK)
  wac = np.sqrt((ka*cia)**2/(1+(ka*cia/wpi)**2))
  wah = np.sqrt( (wpi**2) * (ka*ka * dli*dli * (Te/Ti))/(1+(ka*ka * dli*dli * (Te/Ti))) )
  wl = np.sqrt( (wpe**2) * (1+me/mi) * (1+(3*ka*ka*dl*dl)/(1+me/mi)) )

  kadl = ka*dl


  #
  # coeff1 = 1
  # coeff2 = -2*kadl*(vb/dl)
  # coeff3 = ( (kadl*kadl) * ((vb/dl)*(vb/dl)) ) + ((kadl*kadl)/(1+(kadl*kadl)))*(wpi*wpi)
  #
  #
  # roots = []
  # for i in range(1,len(kadl)):
  #     coeffs = [coeff1, coeff2[i], coeff3[i]]
  #     # coeffs = [coeff1[i], coeff2[i], coeff3[i]]
  #     root = np.roots(coeffs)
  #     roots.append(root)
  # roots = np.array(roots)
  # print(roots[:,0])
  # print(roots[:,1])
  #
  #
  # wbf = np.real(roots[:,0])/wpi
  # wbs = np.real(roots[:,1])/wpi
  #omega_pe

  if norm == "omega_pi":
    omega /= wpi
    Omega /= wpi
  else:
    omega /=wpe
    Omega /=wpe

  wl /= wpe
  # print(wl)
  wac /= wpi
  wah /= wpi
  K *= dl
  # print(wb)

  Z = np.log(np.abs(F))
  # Z = np.abs(F)

  # ==== Figure =============

  ##### FIG SIZE CALC ############
  figsize = np.array([150,150/1.618]) #Figure size in mm
  dpi = 300                         #Print resolution
  ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

  mp.rc('text', usetex=True)
  mp.rc('font', family='sans-serif', size=14, serif='Computer Modern Roman')
  mp.rc('axes', titlesize=14)
  mp.rc('axes', labelsize=14)
  mp.rc('xtick', labelsize=14)
  mp.rc('ytick', labelsize=14)
  mp.rc('legend', fontsize=14)

  fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
  oRange = len(K[:,0]) #for full omega len(K[:,0])
  # oRange = int(oRange/50)
  # print(K[:oRange,:].shape,Omega[:oRange,:].shape,Z[:oRange,:].shape)
  # print(oRange)
  if norm == "omega_pi":
    # oRange = int(oRange/200) #for periodic system in x
    # oRange = int(oRange/10) #for bounded system in x
    plt.pcolor(K[:oRange,:], Omega[:oRange,:], Z[:oRange,:],shading='auto',cmap = 'inferno',vmin=np.min(Z[:oRange,:]),vmax=np.max(Z[:oRange,:])) #np.min(Z[:oRange,:])
    plt.colorbar()
  else:
    # oRange = int(oRange/50)
    plt.pcolor(K[:oRange,:], Omega[:oRange,:], Z[:oRange,:],shading='auto',cmap = 'inferno',vmin=np.min(Z[:oRange,:]),vmax=np.max(Z[:oRange,:]))
    # plt.pcolor(Z[:oRange,:],shading='auto')
    #plt.pcolor(K, Omega, Z,shading='auto',vmin=np.min(Z),vmax=np.max(Z))
    #plt.imshow(K, Omega, Z)
    plt.colorbar()

  if norm == "omega_pi":
    # plt.plot(kadl[1:], wbf, '--w', label="Fast Beam Mode")
    # plt.plot(kadl[1:], wbs, '--k', label="Slow Beam Mode")
    plt.plot(kadl, wac, '--b',label="Acoustic Mode")
    # plt.plot(ka, wb, '--w',label="Beam driven waves")
    leg = ax.legend()
    ax.set_xlabel('$k \lambda_{D}$')
    ax.set_ylabel('$\omega/\omega_{pi}$')
  else:
      plt.axhline(y=1.0, color='b', linestyle='--',label='$\omega_{pe}$')
      plt.plot(kadl, wl, '--w', label="langmuir wave")
      leg = ax.legend()
      ax.set_xlabel('$k \lambda_{D}$')
      ax.set_ylabel('$\omega/\omega_{pe}$')

  ax.set_ylim([0, 2])
  plt.savefig(pjoin(savedir, norm+'_disprel.png'))
  plt.show()
