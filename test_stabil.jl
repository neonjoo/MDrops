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


dir = "pushing_to_limit_langfix_extended"
sourcedir = "/home/andris/sim_data/$dir"
sourcefile = "data00012.jld2"
@load "$sourcedir/$sourcefile" data
points0 = data[1]
faces = data[2]
points0 = Array{Float64}(points0)
faces = Array{Int64}(faces)
t0 = data[3]
sourcefile = "data00013.jld2"
@load "$sourcedir/$sourcefile" data
points1 = data[1]
faces = data[2]
points1 = Array{Float64}(points1)
t1 = data[3]
sourcefile = "data00014.jld2"
@load "$sourcedir/$sourcefile" data
points2 = data[1]
faces = data[2]
points2 = Array{Float64}(points2)
t2 = data[3]

velocities0 = (points1-points0)/(t1-t0)
velocities1 = (points2-points1)/(t2-t1)

using PyPlot
pygui()

fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")

(x, y, z) = [points0[i,:] for i in 1:3]
(vx, vy, vz) = [velocities0[i,:] for i in 1:3]
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


fig2 = figure(figsize=(7,7))
ax2 = fig2[:gca](projection="3d")

(x, y, z) = [points1[i,:] for i in 1:3]
(vx, vy, vz) = [velocities1[i,:] for i in 1:3]
ax2[:scatter](x,y,z, s=2,color="k")

ax2[:quiver](x,y,z,vx,vy,vz, length=0.0001, arrow_length_ratio=0.5)
ax2[:set_title]("Velocities")

ax2[:set_xlim](-2,2)
ax2[:set_ylim](-2,2)
ax2[:set_zlim](-2,2)
ax2[:set_xlabel]("x axis")
ax2[:set_ylabel]("y axis")
ax2[:set_zlabel]("z axis")
fig[:show]()
