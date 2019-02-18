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

points_csv= CSV.read("./meshes/points_sphere.csv", header=0)
faces_csv = CSV.read("./meshes/faces_sphere.csv", header=0)

println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
points1 = copy(points)
points2 = copy(points)
faces = Array{Int64}(faces')
faces2 = copy(faces)
edges = make_edges(faces)
connectivity = make_connectivity(edges)
# H0 = [0, 0, 10]
# eta = 1
# mu = 30
# gamma = 6.9 * 10^-1
# w = 2*pi/50
#
# dt = 0.01
# steps = 50
points = points #.* (4.9 * 10^-1)
#points2 = copy(points)

H0 = [0, 0, 10]


#H0 = 33 .* [0, 0, 1]
#H0 = [0,0,0]
mu = 30
eta = 1
gamma = 6.9 * 10^-1

continue_sim = false
last_step = 0

w = 2*pi/50
t = 0
steps = 10


normals = Normals(points, faces)
for iter in 1:steps
    println("time step $(iter)")

    global points, points1, points2, faces, faces2, normals, edges, connectivity
    #global points2


    (normals, CDE) = make_normals_spline(points, connectivity, edges, normals)
    psi = PotentialSimple(points, faces, mu, H0; normals = normals)
    #psi2 = PotentialSimple(points2, faces, mu, H0; normals = normals2)
    Ht = HtField(points, faces, psi, normals)
    #Ht2 = HtField(points2, faces, psi2, normals2)
    Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
    #Hn_norms2 = NormalFieldCurrent(points2, faces, Ht2, mu, H0; normals = normals2)
    Hn = normals .* Hn_norms'
    #Hn2 = normals2 .* Hn_norms2'

    mup = mu
    # magnitudes squared of the normal force
    Hn_2 = sum(Hn.^2, dims=1)
    #Hn2_2 = sum(Hn2.^2, dims=1)
    # magnitudes squared of the tangential force
    Ht_2 = sum(Ht.^2, dims=1)
    #Ht2_2 = sum(Ht2.^2, dims=1)

    tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2
    #tensorn2 = mup*(mup-1)/8/pi * Hn2_2 + (mup-1)/8/pi * Ht2_2

    # make the force normal to surface
    #tensorn = normals .* tensorn

    velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)
    #velocitiesn_norms2 = InterfaceSpeedZinchenko(points2, faces, tensorn2, eta, gamma, normals2)

    velocitiesn = normals .* velocitiesn_norms'
    #velocitiesn2 = normals2 .* velocitiesn_norms2'


    velocities = velocitiesn
    #velocities2 = make_Vvecs_conjgrad(normals,faces, points, velocitiesn, 1e-6, 120)

    dt = 0.5*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
    #dt2 = 0.4*minl2/maxv2
    println("dt = $(dt)")
    #println("dt2 = $(dt2)")

    points += velocities * dt
    #points2 += velocities2 * dt

end
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)
points1 = active_stabilize(points,faces,CDE,connectivity,normals;maxiters=100)
# for i in 1:size(points,2)
#     for j in 1:size(points,2)
#         edges = make_edges(faces2)
#         connectivity = make_connectivity(edges)
#         println(i," ",j)
#         flip_connectivity!(faces2,i,j,points,connectivity)
#     end
# end
# edges = make_edges(faces2)
# connectivity = make_connectivity(edges)
# (normals2, CDE2) = make_normals_spline(points1, connectivity, edges, normals)
# points2 = active_stabilize(points1,faces2,CDE2,connectivity,normals2;maxiters=100)

# scene = Makie.mesh(points0', faces',color = :white, shading = false,visible = false)
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 1)
# scene = Makie.mesh!(points2', faces',color = :gray, shading = false,visible = true)
# Makie.wireframe!(scene[end][1], color = :blue, linewidth = 1,visible = true)
