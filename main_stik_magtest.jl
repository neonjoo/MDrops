cd("/home/andris/MDrops/")

using JLD2
using StatsBase
using LinearAlgebra
using FastGaussQuadrature
using Optim

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")

## loading the mesh

datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30/"

files = readdir(datadir)
file = files[2+10000]
#println(file)
@load "$datadir/$file" data

points, faces = data[1], data[2]
faces = Array{Int64,2}(faces)
mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                StatsBase.mean(points[2,:]),
                StatsBase.mean(points[3,:]))
points = points .- [mean_x, mean_y, mean_z]

a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)

## setting simulation parameters
H0 = [0., 0., 1.]
mu = 10.
lambda = 10.
Bm = 0.

## calculating the magnetic field on the surface
psi = PotentialSimple(points, faces, normals, mu, H0)
Ht_vec = HtField(points, faces, psi, normals) # a vector
Ht = sqrt.(sum(Ht_vec.^2,dims=1))'
Hn = NormalFieldCurrent(points, faces, normals, Ht_vec, mu, H0) # a scalar

dHn = make_deltaH_normal(points, faces, normals, mu, H0; gaussorder=3)
dHn = make_deltaH_normal(points, faces, normals, mu, H0; gaussorder=3)
Ht_aa = make_H_tangential(points, faces, normals, CDE, H0, dHn; gaussorder = 3)'
Hn_aa = make_H_normal(dHn,mu)'
Ht_aa_flat = make_H_tangential(points, faces, normals, CDE, H0, dHn; gaussorder = 3)'

## testing against theoretical values
using QuadGK
function demag_coefs(a, b, c)

    upper_limit = 2000
    Rq2(q) = (a^2+q) * (b^2+q) * (c^2+q)

    Nx = a*b*c/2 * quadgk(s -> 1/(a^2+s) / sqrt(Rq2(s)), 0, upper_limit)[1]
    Ny = a*b*c/2 * quadgk(s -> 1/(b^2+s) / sqrt(Rq2(s)), 0, upper_limit)[1]
    Nz = a*b*c/2 * quadgk(s -> 1/(c^2+s) / sqrt(Rq2(s)), 0, upper_limit)[1]

    return [Nx, Ny, Nz]
end

function field_theor(a, b, c, mu, H0)

    Ns = demag_coefs(a,b,c)

    Hx = H0[1] / (1 + Ns[1] * (mu-1))
    Hy = H0[2] / (1 + Ns[2] * (mu-1))
    Hz = H0[3] / (1 + Ns[3] * (mu-1))

    return [Hx, Hy, Hz]
end

Hteor = field_theor(a, b, c, mu, H0)
Hn_teor = zeros(size(points,2))
Ht_teor = zeros(size(points,2))
tangs = zeros(size(normals))
for ykey in 1:size(points,2)
    Hn_teor[ykey] = dot(Hteor,normals[:,ykey])
    Ht_teor[ykey] = norm(cross(Hteor,normals[:,ykey]))
end

## maggic
dHn_teor = (mu-1)*Hn_teor'
Ht_aa_teor = make_H_tangential(points, faces, normals, CDE, H0, dHn_teor; gaussorder = 3)'
Ht_aa_teor_flat = make_H_tangential(points, faces, normals, CDE, H0, dHn_teor; gaussorder = 3)'
## plot the results
using Plots
Plots.plot(Hn_teor ./ Hn_teor .- 1,label = "teor",linewidth=2,legend=:bottom)
Plots.plot!(Hn ./ Hn_teor .- 1,label = "old")
Plots.plot!(Hn_aa ./ Hn_teor .- 1,label = "new")

Plots.plot(Hn_teor,label = "teor",linewidth=2,legend=:bottom)
Plots.plot!(Hn,label = "old")
Plots.plot!(Hn_aa,label = "new")

Plots.plot(Ht_teor,label = "teor",linewidth=2,legend=:bottom)
Plots.plot!(Ht,label = "old")
Plots.plot!(Ht_aa_flat,label = "new_flat")
Plots.plot!(Ht_aa,label = "new_curv")
Plots.plot!(Ht_aa_teor_flat,label = "new_flat_teor_Hn")
Plots.plot!(Ht_aa_teor,label = "new_curv_teor_Hn")

Plots.plot(Ht_teor ./ Ht_teor .- 1,label = "teor",linewidth=2,legend=:top)
Plots.plot!(Ht ./ Ht_teor .- 1,label = "old")
Plots.plot!(Ht_aa_flat ./ Ht_teor .- 1,label = "new_flat")
Plots.plot!(Ht_aa ./ Ht_teor .- 1,label = "new_curv")
Plots.plot!(Ht_aa_teor_flat ./ Ht_teor .- 1,label = "new_flat_teor_Hn")
Plots.plot!(Ht_aa_teor ./ Ht_teor .- 1,label = "new_curv_teor_Hn")
