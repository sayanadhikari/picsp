import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
from pylab import *
from spectrum import *



# numerical data file
filename="../output/phi.dat"
time,phi = np.loadtxt(filename, unpack=True)
# disp(data_act.shape)
dt = time[2]-time[1]
sig = phi


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
Omega_pe = np.sqrt((n_0*q*q)/(M_e*eps_0))
Omega_pi = np.sqrt((n_0*q*q)/(M_i*eps_0))
DL       = np.sqrt((TeSI*KB*eps_0)/(n_0*q*q))

print(Omega_pe)

figure(1)
plt.plot(time,phi, linewidth=1,color='blue')
plt.xlabel('Time (s)')
plt.ylabel('$\\phi$ (V)')
plt.title('Potential')

figure(2)
p1 = Periodogram(sig, sampling=1/dt,window='hann',scale_by_freq=False)
p1.run()
p1.plot(norm=True, color='blue', linewidth=0.5)
plt.xscale('symlog')
plt.xlim(1E7,1E11)
plt.xlabel("Frequency / Hz")
plt.ylabel("Power Spectral Density / dB/Hz")
plt.title('Frequency Spectrum')

plt.show()

