using StatsBase
using LinearAlgebra
using FastGaussQuadrature
#using Optim

include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")
#include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")


## making the mesh
points, faces = expand_icosamesh(;R=1,depth=3)
points = Array{Float64}(points)
faces = Array{Int64}(faces)

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)

## setting simulation parameters
H0 = [0., 0., 1.]
mu = 30.
lambda = 10.
Bm = 0.

## calculating the magnetic field on the surface
psi = PotentialSimple(points, faces, normals, mu, H0)
Ht_vec = HtField(points, faces, psi, normals) # a vector
#Ht =
Hn = NormalFieldCurrent(points, faces, normals, Ht_vec, mu, H0) # a scalar
