#from scipy.fft import fft, fftfreq
from scipy.fftpack import fft, fft2
import numpy as np
import matplotlib.pyplot as plt
import plasmapy.dispersion.dispersionfunction
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from PIL import Image
import math
import h5py
from os.path import join as pjoin



EPS = 8.85418782E-12
K = 1.38065E-23
chargeE = 1.602176565E-19
pi = 3.14159265359
massI = 1.673e-27
massE = 9.109e-31

density = 1e12
omega_pe = math.sqrt((chargeE*chargeE*density)/(massE*EPS))
Lambda_D = math.sqrt(EPS*K*11600/(density*chargeE*chargeE))
dt = 0.02 #0.01*(1/omega_pe)
dx = 0.03 #0.01*Lambda_D


path = pjoin("..","output")
print(path)

file_name = 'data'

h5 = h5py.File(pjoin(path,file_name+'.h5'),'r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
# Lz = h5.attrs["Lz"]

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])


Efx = []
for j in range(0,1000,10):
    data_efx = h5["efx/%d"%j]
    efx = data_efx[:,1]
    Efx.append(efx)

Efx = np.array(Efx)

Nx = 64
Ny = 2
Lx = Nx*dx
Omega_l = 0
Omega_h = 150
k_l = 0
k_h = Nx
# x = np.linspace(0,0.784,65)
# y = np.linspace(0,1,100)
# X,Y = np.meshgrid(x, y)
# r,c = (100,65)
# arr = [np.linspace(0,1,100)]*65
# arr = np.ravel(arr)
# arr = np.reshape(arr,(100,65))


Efxf = fft2(Efx)
Exxt = []

F = plasmapy.dispersion.dispersionfunction.plasma_dispersion_func(Efx)
Fd = plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_deriv(Efx)

# def plot_complex(X, Y, Z, N=50):
#     fig, (real_axis, imag_axis) = plt.subplots(1, 2)
#     real_axis.contourf(X, Y, F.real, N)
#     imag_axis.contourf(X, Y, F.imag, N)
#     real_axis.set_title("Real values")
#     imag_axis.set_title("Imaginary values")
#     fig.tight_layout()

# print(arr.shape)
# plot_complex(X,Y,F)
# plt.contourf(F.imag,25)
plt.contourf(Fd[Omega_l:Omega_h,k_l:k_h].real,50)
plt.colorbar()
#
plt.show()
