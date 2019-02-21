using LinearAlgebra
using CSV
using Makie
using StatsBase
using Optim
#using JLD2
#using SurfaceGeometry
#using ElTopo
#using PyPlot
#include("./SurfaceGeometry/dt20L/src/Iterators.jl")
#include("./SurfaceGeometry/dt20L/src/ComplexDS.jl")
include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./sandbox_lang.jl")

points_csv= CSV.read("./meshes/points_bad_ball.csv", header=0)
faces_csv = CSV.read("./meshes/faces_bad_ball.csv", header=0)

println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')

points2 = copy(points)
faces2 = copy(faces)

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
normals, CDE = make_normals_spline(points, connectivity, edges, normals)
closefaces = make_closefaces(faces)

#points2 = active_stabilize(points,faces,CDE,connectivity,normals,deltakoef=0.1,critSc = 0.999,critCdelta = 1.001)#,critSc = 0.75,critCdelta = 1.15)
for i = 1:3
    global faces2,connectivity,points2, CDE, normals
    println("starting stabilization")
    flip_edges!(faces2,connectivity,points2)
    edges = make_edges(faces2)
    normals, CDE = make_normals_spline(points2, connectivity, edges, normals)
    points2 = active_stabilize(points2,faces2,CDE,connectivity,normals,deltakoef=0.1,critSc = 0.999,critCdelta = 1.001)#,critSc = 0.75,critCdelta = 1.15)
end

# scene = Makie.mesh(points', faces',color = :white, shading = false,visible = false)
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 2)
# scene = Makie.mesh(points', faces2',color = :gray, shading = false,visible = true)
# Makie.wireframe!(scene[end][1], color = :blue, linewidth = 2,visible = true)
