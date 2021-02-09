import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
from scipy.spatial import Delaunay, delaunay_plot_2d
from numpy.linalg import norm

folder = "/home/laigars/sim_data/star_5/"
filenames = os.listdir(folder)
filenames.sort()
filename = filenames[-1201] # plot last file

edges = "/home/laigars/mdrops/outer2.jld2"

fig = plt.figure(figsize=(10,5))

for i, file in enumerate([filename]):#filenames[-4:-3]):
    print("working on ",file)    
    f = h5py.File(folder+file, "r")
    
    data = f["data"]
    points = f[data[0]]
    faces = f[data[1]]
    faces = faces[:] -1 # shifting index, because python
    
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    H = f[data[3]][()]
    t = f[data[2]][()]
        
    f = h5py.File(edges,"r")
    outer = f["outer_edges"][:] - 1
    
    ax = fig.add_subplot(1,2,1)
    
    ax.triplot(x, z, faces)    
    ax.scatter(x, z, marker=".", c="black", zorder=2)
    
    ax2 = fig.add_subplot(1,2,2)
    
    ax2.triplot(x, z, faces)    
    ax2.scatter(x, z, marker=".", c="black", zorder=2)
    ax2.set_xlim([-1.9, -1.5])
    ax2.set_ylim([-3.3, -2.5])
    
    for nodes in outer:
        n1, n2 = nodes
        xs, zs = x[n1], z[n1]
        xe, ze = x[n2], z[n2]
        #print(norm(xe-xs, ))
        ax.plot([xs, xe], [zs, ze], "r", lw=2)
        ax2.plot([xs, xe], [zs, ze], "r", lw=2)

        
        
        
        
        
        
        