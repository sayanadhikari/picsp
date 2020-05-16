import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec


# numerical data file
filename="potential.dat"

data_act = np.loadtxt(filename, delimiter="\t", skiprows=2, usecols=[0,2])

data1 = data_act[:,0] #time data
data2 = data_act[:,1] #potential data


t = data1
dt = data1[2]-data1[1]
sig = data2

plt.subplot(211)
plt.plot(t, sig)
plt.subplot(212)
plt.psd(sig, 512, 1 / dt)
plt.xscale('symlog')

plt.show()
