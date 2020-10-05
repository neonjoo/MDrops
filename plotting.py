import h5py
from matplotlib import pyplot as plt
import matplotlib.tri as mtri
import os
import numpy as np

# read file
folder = "/home/laigars/sim_data/star_4/"
filenames = os.listdir(folder)
filenames.sort()
filename = filenames[-2] # plot last file
#framenumber = 1
#filename = filenames[framenumber+1] # plot framenumber file
print("showing ",filename)

f = h5py.File(folder+filename, "r")
#f = h5py.File("/home/andris/Desktop/before_dis.jld2", "r")

data = f["data"]
points = f[data[0]]
faces = f[data[1]]
faces = faces[:] -1 # shifting index, because python

x = points[:,0]
y = points[:,1]
z = points[:,2]

t = f[data[2]]#[()]

# Plotting


fig = plt.figure()
ax = fig.gca(projection='3d')
triang = mtri.Triangulation(x, y, faces)
ax.plot_trisurf(triang, z,lw=0.2, edgecolor="black", color="white")

axlim = 3
ax.set_xlim3d([-axlim, axlim])
ax.set_ylim3d([-axlim, axlim])
ax.set_zlim3d([-axlim, axlim])
ax.view_init(0,0)
plt.show()
