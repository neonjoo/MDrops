using Pkg
pkg"activate ."
pkg"resolve"

using JLD2
using StatsBase
using LinearAlgebra
using FastGaussQuadrature
using Optim

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")

## making the mesh
points, faces = expand_icosamesh(;R=1,depth=2)
points = Array{Float64}(points)
faces = Array{Int64}(faces)




edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
#(normals, CDE, AB) = make_normals_parab(points, connectivity, edges, normals)

# #
# ## setting simulation parameters
H0 = [0., 0., 1.]
mu = 10.
lambda = 10.
Bm = 3.

# calculating the magnetic field on the surface
@time psi = PotentialSimple(points, faces, normals, mu, H0)
@time psi_par = PotentialSimple_par(points, faces, normals, mu, H0)

Ht_vec = HtField_par(points, faces, psi, normals) # a vector
# Ht = sqrt.(sum(Ht_vec.^2,dims=1))'
@time Hn_norms = NormalFieldCurrent_par(points, faces, normals, Ht_vec, mu, H0) # a scalar
Hn = normals .* Hn_norms'

Hn_2 = sum(Hn.^2, dims=1)
Ht_2 = sum(Ht_vec.^2, dims=1)

@time v_phys = make_magvelocities_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)

vvecs = v_phys
#
# println(333333)
# @time v_stab = make_Vvecs_conjgrad(normals, faces, points, vvecs, 1e-6, 500)
# @time v_stab_par = make_Vvecs_conjgrad_par(normals, faces, points, vvecs, 1e-6, 500)
# @time v_stab_par2 = make_Vvecs_conjgrad_par2(normals, faces, points, vvecs, 1e-6, 500)
