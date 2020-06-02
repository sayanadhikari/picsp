import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pylab import *
# numerical data file
filename1="../output/phase_space/i20000.dat"
filename2="../output/phase_space/e20000.dat"

x1,v1 = np.loadtxt(filename1, unpack=True)
x2,v2 = np.loadtxt(filename2, unpack=True)
# disp(data_act.shape)

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


figure(1)
plt.scatter(x1/DL,v1,s=1,marker='.')
plt.xlabel("$x/\\lambda_D$")
plt.ylabel("$V_x$")
plt.title('Ion Phase Space')
plt.tight_layout()


figure(2)
plt.scatter(x2/DL,v2,s=1,marker='.')
plt.xlabel("$x/\\lambda_D$")
plt.ylabel("$V_x$")
plt.title('Electron Phase Space')
plt.tight_layout()

plt.show()

