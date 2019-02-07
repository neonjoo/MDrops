using LinearAlgebra
using CSV
using JLD2
using Plots
include("./SurfaceGeometry/dt20L/src/Iterators.jl")
#using SurfaceGeometry

using PyPlot
#
pygui(true)
pyplot()

datadir = "./data/stikuts1"

fig = figure()#figsize=(7,7))
#ax = fig[:gca](projection="3d")
#@load "./data/simul/points50.jld2" points2
#allp = points2[:,:]

for i in [2, 9000]

    ax1=fig[:gca](projection="3d")#, aspect="equal")
    #ax1[:plot_wireframe](x,y,z)
    #ax1[:set_axis_off]()
    #ax2=fig[:add_subplot](2,1,2,projection="3d", aspect="equal")
    #ax2[:plot_wireframe](x,y,z)
    #ax2[:set_axis_off]()
    println(i)
    # if i % 500 == 0
    global data
    @load "$datadir/data$(lpad(i,5,"0")).jld2" data
    #global allp
    #allp = cat(dims=3, allp, points2)
    #global x,y,z
    points2, H0 = data[1], data[3]
    (x, y, z) = [points2[j,:] for j in 1:3]
    # end
    ax1[:scatter](x,y,z, s=3, label="step $i, H=$(round.(H0, digits=3))")#,color="k")
    lowlim, uplim = (-5, 5) .* 10^-4
    ax1[:set_xlim](lowlim, uplim)
    ax1[:set_ylim](lowlim, uplim)
    ax1[:set_zlim](lowlim, uplim)
    ax1[:set_xlabel]("x axis")
    ax1[:set_ylabel]("y axis")
    ax1[:set_zlabel]("z axis")
    ax1[:legend]()
    # @load "./data/simul/points$i.jld2" points2
    #
    # (x, y, z) = [points2[j,:] for j in 1:3]
    # # end
    # ax2[:scatter](x,y,z, s=3, label="not-stabil")#,color="k")
    # ax2[:set_xlim](-2,2)
    # ax2[:set_ylim](-2,2)
    # ax2[:set_zlim](-2,2)
    # ax2[:set_xlabel]("x axis")
    # ax2[:set_ylabel]("y axis")
    # ax2[:set_zlabel]("z axis")
    # ax2[:legend]()
end

fig[:show]()

#ax[:quiver](x,y,z,vnx,vny,vnz, length=2, arrow_length_ratio=0.5)

#ax[:quiver](x,y,z,Hnx,Hny,Hnz, length=0.3, arrow_length_ratio=0.5)
#ax[:quiver](x,y,z,Htx,Hty,Htz, length=0.3, arrow_length_ratio=0.5)
#ax[:quiver](x,y,z,Hx,Hy,Hz, length=0.3)
# #
# ax[:set_xlim](-2,2)
# ax[:set_ylim](-2,2)
# ax[:set_zlim](-2,2)
# ax[:set_xlabel]("x axis")
# ax[:set_ylabel]("y axis")
# ax[:set_zlabel]("z axis")
# ax[:legend]()
# fig[:show]()

#
# (Htx, Hty, Htz) = Ht[1,:], Ht[2,:], Ht[3,:]
# #velocitiesn = all_normals .* velocitiesn_norm'
# (vnx, vny, vnz) = [velocitiesn[i,:] for i in 1:3]
# (Hnx, Hny, Hnz) = [Hn[i,:] for i in 1:3]
# Hn_teor = normals .* Hn_teor'
# #Ht_teor = normals .* Ht_teor'
# (Hnx_teor, Hny_teor, Hnz_teor) = [Hn_teor[i,:] for i in 1:3]
# (Htx_teor, Hty_teor, Htz_teor) = [Ht_teor[i,:] for i in 1:3]
# (fx, fy, fz) = [tensorn[i,:] for i in 1:3]
# (nx, ny, nz) = [normals[i,:] for i in 1:3]
# (Hx, Hy, Hz) = (Htx+Hnx, Hty+Hny, Htz+Hnz)
#
# fig = figure(figsize=(7,7))
# ax = fig[:gca](projection="3d")
#
# (x, y, z) = [points2[i,:] for i in 1:3]
# ax[:scatter](x,y,z, s=2,color="k")
# #(xp, yp, zp) = points[:,418]
# #ax[:scatter](xp, yp, zp, s=10,color="r")
#
# # for v in NeighborVertices(418, faces)
# #     (xp, yp, zp) = points[:,v]
# #     ax[:scatter](xp, yp, zp, s=10,color="g")
# # end
#
# #ax[:set_title]("Normals")
# #ax[:quiver](x,y,z,nx,ny,nz, length=0.2, arrow_length_ratio=0.5)
#
# ax[:quiver](x,y,z,fx,fy,fz, length=5, arrow_length_ratio=0.5)
# ax[:set_title]("Normal forces")
#
# #ax[:quiver](x,y,z,vnx,vny,vnz, length=5, arrow_length_ratio=0.5)
# #ax[:set_title]("Normal velocities")
#
# #ax[:quiver](x,y,z,Hnx,Hny,Hnz, length=0.3, arrow_length_ratio=0.5, color="blue")
# #ax[:quiver](x,y,z,Hnx_teor,Hny_teor,Hnz_teor, length=0.4, arrow_length_ratio=0.5,color="blue")
# #ax[:set_title]("Normal H component, external H // z")
#
# #ax[:quiver](x,y,z,Htx,Hty,Htz, length=0.3, arrow_length_ratio=0.5, color="blue")
# #ax[:quiver](x,y,z,Htx_teor,Hty_teor,Htz_teor, length=0.3, arrow_length_ratio=0.5,color="red")
# #ax[:set_title]("Tangential H")
#
# #ax[:quiver](x,y,z,Hx,Hy,Hz, length=0.3)
# ax[:set_xlim](-2,2)
# ax[:set_ylim](-2,2)
# ax[:set_zlim](-2,2)
# ax[:set_xlabel]("x axis")
# ax[:set_ylabel]("y axis")
# ax[:set_zlabel]("z axis")
# fig[:show]()
#
