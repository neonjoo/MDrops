using LinearAlgebra
using CSV
using Makie
using StatsBase
using Optim
using JLD2
#using SurfaceGeometry
#using ElTopo
#using PyPlot
#include("./SurfaceGeometry/dt20L/src/Iterators.jl")
#include("./SurfaceGeometry/dt20L/src/ComplexDS.jl")
include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry
#include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./sandbox_lang.jl")
include("./physics_functions.jl")


dir = "2019-03-12/7"
sourcedir = "/home/andris/sim_data/$dir"
fileNo = 1
sourcefile = "data$(lpad(fileNo,5,"0")).jld2"
@load "$sourcedir/$sourcefile" data
points = data[1]
faces = data[2]
points = Array{Float64}(points)
faces = Array{Int64}(faces)
t = data[3]
sourcefile = "data$(lpad(fileNo+1,5,"0")).jld2"
@load "$sourcedir/$sourcefile" data
points1 = data[1]
faces1 = data[2]
points1 = Array{Float64}(points1)
faces1 = Array{Int64}(faces1)
t1 = data[3]

velocities = (points1-points)/(t1-t)

using PyPlot
pygui()

fig = PyPlot.figure(figsize=(7,7))
#ax = fig[:gca](projection="3d")
#ax = getproperty(fig,:gca)
#ax(projection = "3d")

(x, y, z) = [points[i,:] for i in 1:3]
(vx, vy, vz) = [velocities[i,:] for i in 1:3]
PyPlot.plot3D(x,y,z;color="k")
ax = fig[:gca]
ax[:scatter](x,y,z, s=2,color="k")

ax[:quiver](x,y,z,vx,vy,vz, length=0.0001, arrow_length_ratio=0.5)
ax[:set_title]("Velocities")

ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()
