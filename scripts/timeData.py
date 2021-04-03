#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse
import matplotlib as mp
import os
# import plotly.graph_objects as go
#========= Configuration ===========
parser = argparse.ArgumentParser(description='Grid Data Animator PICSP')
parser.add_argument('-p', default="energy", type=str, help='Name of the parameter (energy)')

args = parser.parse_args()

param = args.p

DIR ="../output/"

file_name = "data"#"rhoNeutral" #"P"


#========== Figure Directory Setup =============
figPath  = "figures"  # DO NOT CHANGE THE PATH
if os.path.exists(figPath):
    print("figure directory exists. Existing figures will be replaced.")
else:
    os.mkdir(figPath)

#=========== Data loading from HDF5 file ================
h5 = h5py.File(DIR+file_name+'.h5','r')

# Attribute load from data h5
dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

# dataset index
data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

# Data load from data h5
data = h5["timedata/"+param]

# ==== Figure =============

##### FIG SIZE CALC ############
figsize = np.array([77,77/1.618]) #Figure size in mm
dpi = 300                         #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
mp.rc('text', usetex=True)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)

ax.plot(data_num,data[:-1,0],label="$KE_{i}$")
ax.plot(data_num,data[:-1,1],label="$KE_{e}$")
leg = ax.legend()
ax.set_xlabel('$TimeSteps\,[\mathrm{m}]$')
if ("energy" in param ):
    ax.set_ylabel('$Energy\,[\mathrm{A.U.}]$')
    figName = "energy"
plt.savefig(figPath+"/"+figName+"-PICSP",dpi=dpi)
plt.show() #After savefig
