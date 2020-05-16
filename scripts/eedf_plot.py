import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pylab import *
# numerical data file
filename="../output/phase_space/i1000.dat"

data1,data2 = np.loadtxt(filename, unpack=True)
# disp(data_act.shape)


#velocity data within a range
range1 = 0.001
range2 = 0.002



#Empty data array to store velocity data for a location range
data = np.empty([len(data1), 1], dtype=float)
data[:] = np.NaN

#storing data for the range given
for i in range(len(data1)):
    if (data1[i] >= range1) and (data1[i] <= range2):
        data[i]= data2[i]

#Conversion of velocities to energy
#data in eV
mass = 9.1E-31 #electron
charge_electron = 1.6e-19
KE_in_EV_coeff = 0.5*mass/charge_electron
datasqr = data*data
data_EV = KE_in_EV_coeff*datasqr

figure(1)
sns.distplot(data_EV[:,0], hist=False, kde=True, color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'shade': True, 'linewidth': 2})
plt.xlabel('Energy (eV)')
plt.ylabel('A.U.')
plt.title('Electron Energy Distribution Function (EEDF)')


figure(2)
sns.distplot(data[:,0], hist=False, kde=True, color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'shade': True, 'linewidth': 2})
plt.xlabel('Velocity (m/s)')
plt.ylabel('A.U.')
plt.title('Electron Velocity Distribution Function (EVDF)')
plt.show()
