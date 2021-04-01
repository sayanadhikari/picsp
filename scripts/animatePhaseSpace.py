#!/usr/bin/ python

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
# import plotly.graph_objects as go
#========= Configuration ===========
show_anim = True
save_anim = False
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
        ax1.set_title('Electron Phase Space (TimeSteps = %d'%i+')')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$v_x$")
        ax1.set_xlim([0, Lx])

        ax2.cla()
        img2 = ax2.scatter(dataix,dataivx,marker='.',color='r',alpha=1.0,s=10)
        ax2.set_title('Ion Phase Space (TimeSteps = %d'%i+')')
        ax2.set_xlabel("$y$")
        ax2.set_ylabel("$v_y$")
        ax2.set_xlim([0, Ly])
        # ax1.set_ylim([-1, 1])
        # ax1.set_zlim([-Lz, Lz])



if (show_anim == True):
    fig,(ax1,ax2) = plt.subplots(2,1, figsize=(8, 8))
    # fig = plt.figure(figsize=(6, 6))
    # ax1 = plt.axes(projection ="3d")
    ani = animation.FuncAnimation(fig,animate,frames=len(data_num),interval=interval*1e+3,blit=False)
    # ani.save('phase_space.gif',writer='imagemagick')
    plt.show()
    if(save_anim == True):
        try:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=(1/interval), metadata=dict(artist='Me'), bitrate=1800)
        except RuntimeError:
            print("ffmpeg not available trying ImageMagickWriter")
            writer = animation.ImageMagickWriter(fps=(1/interval))
        ani.save('animation2d.mp4')
