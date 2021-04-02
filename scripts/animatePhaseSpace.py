#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
import argparse
# import plotly.graph_objects as go
#========= Configuration ===========


parser = argparse.ArgumentParser(description='Phase - Space Data Animator PICSP')
parser.add_argument('-a', default=True, type=bool, help='Show Animation (True/False)')
parser.add_argument('-s', default=False, type=bool, help='Save Animation (True/False)')
args = parser.parse_args()

show_anim = args.a
save_anim = args.s

interval = 0.001#in seconds

DIR ="../output/"

file_name = "data"#"rhoNeutral" #"P"

h5 = h5py.File(DIR+file_name+'.h5','r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]
# Lz = h5.attrs["Lz"]

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])


data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

if (show_anim == True):
    def animate(i):
        #======Electron Data=========
        datae = h5["particle.e/%d"%data_num[i]]
        dataex = datae[:,0]
        dataey = datae[:,1]
        dataevx = datae[:,2]
        dataevy = datae[:,3]

        #======Ion Data=========
        datai = h5["particle.i/%d"%data_num[i]]
        dataix = datai[:,0]
        dataiy = datai[:,1]
        dataivx = datai[:,2]
        dataivy = datai[:,3]

        ax1.cla()
        img1 = ax1.scatter(dataex,dataevx,marker='.',color='b',alpha=1.0,s=10)
        ax1.set_title('Electron Phase Space (TimeSteps = %d'%(i*dp)+')')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$v_x$")
        ax1.set_xlim([0, Lx])

        ax2.cla()
        img2 = ax2.scatter(dataix,dataivx,marker='.',color='r',alpha=1.0,s=10)
        ax2.set_title('Ion Phase Space (TimeSteps = %d'%(i*dp)+')')
        ax2.set_xlabel("$y$")
        ax2.set_ylabel("$v_y$")
        ax2.set_xlim([0, Ly])
        # ax1.set_ylim([-1, 1])
        # ax1.set_zlim([-Lz, Lz])




fig,(ax1,ax2) = plt.subplots(2,1, figsize=(8, 6))
ani = animation.FuncAnimation(fig,animate,frames=len(data_num),interval=interval*1e+3,blit=False)
if (show_anim == True):
    plt.show()
if(save_anim == True):
    try:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=(1/interval), metadata=dict(artist='Me'), bitrate=1800)
    except RuntimeError:
        print("ffmpeg not available trying ImageMagickWriter")
        writer = animation.ImageMagickWriter(fps=(1/interval))
    ani.save('PICSP-phasespace.mp4')
