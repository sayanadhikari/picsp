import numpy as np
import h5py
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Grid Data Viewer PICSP')
parser.add_argument('-p', default="efx", type=str, help='Name of the parameter (phi/den.e/den.i)')
parser.add_argument('-t', default=200, type=int, help='timestep to view (e.g. 100)')

args = parser.parse_args()


param = args.p
data_num = args.t

DIR ="../output/"

file_name = "data"#"rhoNeutral" #"P"


h5 = h5py.File(DIR+file_name+'.h5','r')

Lx = h5.attrs["Lx"]
Ly = h5.attrs["Ly"]

Nx = int(h5.attrs["Nx"])
Ny = int(h5.attrs["Ny"])

dp   =  int(h5.attrs["dp"])
Nt   =  int(h5.attrs["Nt"])
den_norm_factor = h5.attrs["den_norm_factor"]

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X, Y = np.meshgrid(x, y)

data = h5[param+"/%d"%data_num]
# data = np.transpose(data)
if ("den" in param ):
    data = data*den_norm_factor
# Creating pandas dataframe from numpy array
# dataset = pd.DataFrame({'Column1': data[:, 0], 'Column2': data[:, 1]})
dataset = pd.DataFrame(data)
print("Parameter: "+param+" Timestep: %d"%data_num)
print(dataset)
