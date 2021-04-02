#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
# import plotly.graph_objects as go
#========= Configuration ===========
parser = argparse.ArgumentParser(description='Grid Data Animator PICSP')
parser.add_argument('-p', default="phi", type=str, help='Name of the parameter (phi/den.e/den.i)')
parser.add_argument('-a', default=True, type=bool, help='Show Animation (True/False)')
parser.add_argument('-s', default=False, type=bool, help='Save Animation (True/False)')
parser.add_argument('-d', default=False, type=bool, help='3D Animation (True/False)')
args = parser.parse_args()

param = args.p
show_anim = args.a
save_anim = args.s
Vis3D = args.d
interval = 0.001#in seconds

DIR ="../output/"

file_name = "data"#"rhoNeutral" #"P"

h5 = h5py.File(DIR+file_name+'.h5','r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]

Nx = int(h5.attrs["Nx"])
Ny = int(h5.attrs["Ny"])

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X, Y = np.meshgrid(x, y)

# dataset index
data_num = np.arange(start=0, stop=Nt, step=dp, dtype=int)

maxPhi = np.max(h5[param+"/%d"%Nt]);
minPhi = np.min(h5[param+"/%d"%Nt]);

if (show_anim == True):
    def animate(i):
        #======Potential Data=========
        data = h5[param+"/%d"%data_num[i]]


        ax1.cla()
        if Vis3D == True:
            img1 = ax1.plot_surface(X,Y,data)
            ax1.set_zlim([minPhi, maxPhi])
        else:
            img1 = ax1.contourf(X,Y,data)
        if ("phi" in param ):
            ax1.set_title('Potential (TimeSteps = %d'%(i*dp)+')')
        elif ("den" in param ):
            if ("i" in param ):
                ax1.set_title('Ion Density (TimeSteps = %d'%(i*dp)+')')
            else:
                ax1.set_title('Electron Density (TimeSteps = %d'%(i*dp)+')')
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("$y$")
        ax1.set_xlim([0, Lx])
        ax1.set_ylim([0, Ly])
        cax.cla()
        fig.colorbar(img1, cax=cax)





if (show_anim == True):
    fig,ax1 = plt.subplots(1,1, figsize=(6, 6))
    div = make_axes_locatable(ax1)
    cax = div.append_axes('right', '5%', '5%')
    data = h5["phi/%d"%data_num[0]]
    if Vis3D == True:
        fig = plt.figure(figsize=(6, 6))
        ax1 = plt.axes(projection ="3d")
        img1 = ax1.plot_surface(X,Y,data)
    else:
        img1 = ax1.contourf(X,Y,data)
    cbar = fig.colorbar(img1,cax=cax)
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
