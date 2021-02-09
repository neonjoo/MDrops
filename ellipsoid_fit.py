import numpy as np
from   numpy.linalg import eig, inv
import h5py
import os
from matplotlib import pyplot as plt
from scipy import optimize

def ellipsoid_fit(x,y,z):
#    x = X[:, 0]
#    y = X[:, 1]
#    z = X[:, 2]
    D = np.array([x * x + y * y - 2 * z * z,
                 x * x + z * z - 2 * y * y,
                 2 * x * y,
                 2 * x * z,
                 2 * y * z,
                 2 * x,
                 2 * y,
                 2 * z,
                 1 - 0 * x])
    d2 = np.array(x * x + y * y + z * z).T # rhs for LLSQ
    u = np.linalg.solve(D.dot(D.T), D.dot(d2))
    a = np.array([u[0] + 1 * u[1] - 1])
    b = np.array([u[0] - 2 * u[1] - 1])
    c = np.array([u[1] - 2 * u[0] - 1])
    v = np.concatenate([a, b, c, u[2:]], axis=0).flatten()
    A = np.array([[v[0], v[3], v[4], v[6]],
                  [v[3], v[1], v[5], v[7]],
                  [v[4], v[5], v[2], v[8]],
                  [v[6], v[7], v[8], v[9]]])

    center = np.linalg.solve(- A[:3, :3], v[6:9])

    translation_matrix = np.eye(4)
    translation_matrix[3, :3] = center.T

    R = translation_matrix.dot(A).dot(translation_matrix.T)

    evals, evecs = np.linalg.eig(R[:3, :3] / -R[3, 3])
    evecs = evecs.T

    radii = np.sqrt(1. / np.abs(evals))
    radii *= np.sign(evals)

    return center, evecs, radii, v

def ellipsoid_plot(center, radii, rotation, ax, plot_axes=False, cage_color='b', cage_alpha=0.2):
    """Plot an ellipsoid"""
        
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    
    # cartesian coordinates that correspond to the spherical angles:
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    # rotate accordingly
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) + center

    if plot_axes:
        # make some purdy axes
        axes = np.array([[radii[0],0.0,0.0],
                         [0.0,radii[1],0.0],
                         [0.0,0.0,radii[2]]])
        # rotate accordingly
        for i in range(len(axes)):
            axes[i] = np.dot(axes[i], rotation)

        # plot axes
        for p in axes:
            X3 = np.linspace(-p[0], p[0], 100) + center[0]
            Y3 = np.linspace(-p[1], p[1], 100) + center[1]
            Z3 = np.linspace(-p[2], p[2], 100) + center[2]
            ax.plot(X3, Y3, Z3, color=cage_color)

    # plot ellipsoid
    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cage_color, alpha=cage_alpha)

    
# datafolder = "/home/andris/sim_data/"
# folders = os.listdir(datafolder)

# cs = []
# As = []
# bs = []
# Bms = []
# ts = []
# startframe = 1
# for folder in folders:
#     if "quasi_static_test2" in folder:
#         filenames = os.listdir(datafolder + folder)
#         filenames.sort()
        
#         for filename in filenames[startframe +1:-2]:
#             print(filename)

#             f = h5py.File(datafolder + folder+"/"+filename, "r")
#             data = f["data"]
            
#             t = f[data[2]][()]
            
#             points = f[data[0]]
#             x = points[:,0]
#             y = points[:,1]
#             z = points[:,2]
            
#             center, evecs, axes, v = ellipsoid_fit(x,y,z)
            
#             a = axes[0]
#             b = axes[1]
#             c = axes[2]
      
#             Bm = f[data[12]][()]
            
#             cs.append(c)
#             As.append(a)
#             bs.append(b)
#             Bms.append(Bm)
#             ts.append(t)

def Bm_vs_K(K, mu):
    bm = np.pi*( np.sqrt(K**2-1)*(-1 + 2*(K**2-1)/(mu-1)) + K**2 * np.arctan(np.sqrt(K**2-1)))**2 * (
            2*np.sqrt(K**2-1) * K*(1+2*K**2) + (1-4*K**2) * np.log( (K+np.sqrt(K**2-1))/(K-np.sqrt(K**2-1)) )) / (
            K**(7/3) * (K**2-1)**2 * ( -3*np.sqrt(K**2-1) + (2+K**2)*np.arctan(np.sqrt(K**2-1)) )
           )
    return bm

def K_vs_c(c):
    # K = a/c
    # a*b*c = 1
    # a=b
    # a^2 * c = 1 -> a = sqrt(1/c) -> K = 1/c^(3/2)
    return 1/c**(3/2)

def K_vs_a_or_b(a):
    # K = a/c
    # a*b*c = 1
    # a=b
    # a^2 * c = 1 -> a = sqrt(1/c) -> K = 1/c^(3/2)
    # a^2 * c = 1 -> c = 1/a^2 -> K = a/(1/a^2) = a^3
    return a**3

# plt.plot(ts,As)
# plt.plot(ts,bs)
# plt.plot(ts,cs)
# plt.show()
#plt.plot(ts,Bms)

def eqlib_inds(Bms):
    inds = []
    Bm0 = 1
    for (ind,Bm) in enumerate(Bms):
        if Bm > Bm0:
            inds.append(ind-1)
            Bm0 = Bm
    
    inds.append(len(Bms)-1)
    return inds


#eqinds = eqlib_inds(Bms)


#%%
teorC = np.linspace(1/4,0.99,100)
teorCK = K_vs_c(teorC)
teorCBm = Bm_vs_K(teorCK, 6)

teorB = np.linspace(1.01,2,100)
teorBK = K_vs_a_or_b(teorB)
teorBBm = Bm_vs_K(teorBK, 6)

#%%
# plt.scatter(np.array(Bms)[eqinds], np.array(As)[eqinds])
# plt.scatter(np.array(Bms)[eqinds], np.array(bs)[eqinds])
# plt.scatter(np.array(Bms)[eqinds], np.array(cs)[eqinds])
#plt.plot(np.array(Bms)[eqinds], np.array(cs)[eqinds] * np.array(As)[eqinds] * np.array(bs)[eqinds])
#plt.plot(np.array(Bms)[eqinds], np.ones(np.array(Bms)[eqinds].shape))
# plt.plot(teortriBms, teortriAs)
# plt.plot(teortriBms, teortriBs)
# plt.plot(teortriBms, teortriCs)

plt.plot(teorBBm, teorB)
plt.plot(teorBBm, teorB)
plt.plot(teorBBm, teorC)

plt.show()

#%%
plt.scatter(np.array(Bms)[eqinds], np.array(As)[eqinds])
plt.scatter(np.array(Bms)[eqinds], np.array(bs)[eqinds])
plt.scatter(np.array(Bms)[eqinds], np.array(cs)[eqinds])
#plt.plot(np.array(Bms)[eqinds], np.array(cs)[eqinds] * np.array(As)[eqinds] * np.array(bs)[eqinds])
#plt.plot(np.array(Bms)[eqinds], np.ones(np.array(Bms)[eqinds].shape))
plt.plot(teortriBms, teortriAs)
plt.plot(teortriBms, teortriBs)
plt.plot(teortriBms, teortriCs)

plt.plot(teorBBm, teorB)
plt.plot(teorBBm, teorB)
plt.plot(teorCBm, teorC)

plt.xlim(0, 30)
plt.ylim(0.5, 1.5)

plt.show()


teortriBms = []
teortriCs = [] 
for line in open("cs_forw_mu6_bm1-100.dat", 'r'):
    item = line.rstrip().split() # strip off newline and any other trailing whitespace
    teortriBms.append(float(item[0]))
    teortriCs.append(float(item[1]))

teorC = np.linspace(0.2,0.99,100)
teorK = K_vs_c(teorC)
teorBm = Bm_vs_K(teorK, 10)


#plt.plot(teorBm,teorC)
#plt.plot(teortriBms,teortriCs)
#plt.scatter(Bms,cs)
#plt.show()

teortriBms = []
teortriAs = [] 
for line in open("as_forw_mu6_bm1-100.dat", 'r'):
    item = line.rstrip().split() # strip off newline and any other trailing whitespace
    teortriBms.append(float(item[0]))
    teortriAs.append(float(item[1]))
    
#plt.plot(teortriBms,teortriAs)
#plt.scatter(Bms,As)
#plt.show()

teortriBms = []
teortriBs = [] 
for line in open("bs_forw_mu6_bm1-100.dat", 'r'):
    item = line.rstrip().split() # strip off newline and any other trailing whitespace
    teortriBms.append(float(item[0]))
    teortriBs.append(float(item[1]))
#    
#plt.plot(teortriBms,teortriBs)
#plt.scatter(Bms,bs)
#plt.show()
#
#print(teortriAs[779])
#print(teortriBs[779])
#print(teortriCs[779])
#print(teortriAs[779]*teortriBs[779]*teortriCs[779])