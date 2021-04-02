#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse
# import plotly.graph_objects as go
#========= Configuration ===========
parser = argparse.ArgumentParser(description='Grid Data Animator PICSP')
parser.add_argument('-p', default="energy", type=str, help='Name of the parameter (energy)')

args = parser.parse_args()

param = args.p

DIR ="../output/"

file_name = "data"#"rhoNeutral" #"P"

h5 = h5py.File(DIR+file_name+'.h5','r')


dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

# dataset index
data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)



#====== Data=========
data = h5["timedata/"+param]

# ==== Figure =============
fig,ax = plt.subplots(1,1, figsize=(6, 4))
ax.plot(data_num,data[:-1,0],label="$KE_{i}$")
ax.plot(data_num,data[:-1,1],label="$KE_{e}$")
plt.legend()
ax.autoscale(enable=True, axis='x', tight=True)
ax.set_xlabel("$TimeSteps$")
if ("energy" in param ):
    ax.set_ylabel("$Energy$")

plt.show()
