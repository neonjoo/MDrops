using Pkg
pkg"activate ."
pkg"resolve"

using Plots
#using Makie
using JLD2
using FileIO
using Optim
using CSV

p = Plots
p.pyplot()

include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")



points_csv= CSV.read("./meshes/points_sphere.csv", header=0)
faces_csv = CSV.read("./meshes/faces_sphere.csv", header=0)
#fields = CSV.read("/home/laigars/sim_data/field.csv", header=0)[1] * 10 # mT -> Oe
#times = CSV.read("/home/laigars/sim_data/time.csv", header=0)[1]

#points_csv= CSV.read("./meshes/points_ellipse_manyN.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_ellipse_manyN.csv", header=0)
println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')


H0 = [0., 0., 1.]
mu = 30

# for mu=30 hist jump @ Bm=3.69
Bm = 3
eta= "gatyt"
R0 = 21.5 * 100/480 * 1e-4 # um to cm for cgs
R0 = 1
#lambda = 10.1
lambda = 7.6
gamma = H0[3]^2 * R0 / Bm
#gamma = 8.2 * 1e-4
#gamma = 7.7 * 1e-4 # from fitted exp data with mu=34

reset_vmax = true
last_step = 0
t = 0
dt = 0.1
steps = 1

#points, faces = data[1], data[2]
normals = Normals(points, faces)


for i in 1:steps
    println("------------------------------------------------------------------------------------------------- Step ($i)$(i+last_step)")
    global points, faces, connectivity, normals, all_vs, velocities, velocities2
    global t, H0

    edges = make_edges(faces)
    connectivity = make_connectivity(edges)
    normals, CDE = make_normals_spline(points, connectivity, edges, normals)

    psi = PotentialSimple(points, faces, mu, H0; normals = normals)
    Ht = HtField(points, faces, psi, normals)
    Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
    Hn = normals .* Hn_norms'

    gamma = H0[3]^2 * R0 / Bm
    println("gamma = $gamma")
    #mup = mu
    # magnitudes squared of the normal force
    Hn_2 = sum(Hn.^2, dims=1)
    # magnitudes squared of the tangential force
    Ht_2 = sum(Ht.^2, dims=1)


    println("Bm = $Bm")
    # "_2" has early exit for lambda=1
    println("velocity:")
    @time velocities = make_magvelocities(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
    @time velocities2 = make_magvelocities_no_Wie(points, normals, lambda, Bm, mu, Hn_2, Ht_2)

    #@time velocities = make_Vvecs_conjgrad(normals,faces, points, velocities, 1e-6, 500)





end

p.plot(velocities[1,:], label="wie")
p.plot!(velocities2[1,:], label="no wie")
