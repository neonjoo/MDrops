using Pkg
pkg"activate ."
pkg"resolve"

using LinearAlgebra
using CSV
using JLD2
#using Makie
using StatsBase
using Optim
#using ElTopo


include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry
#include("./SurfaceGeometry/dt20L/src/StabilisationMethods/stabilisationV2.jl")

#include("./SurfaceGeometry/dt20L/src/Properties.jl")

#include("./SurfaceGeometry/dt20L/src/Utils.jl")

include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")

points_csv= CSV.read("./meshes/points_sphere.csv", header=0)
faces_csv = CSV.read("./meshes/faces_sphere.csv", header=0)

#points_csv= CSV.read("./meshes/points_ellipse_manyN.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_ellipse_manyN.csv", header=0)
println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')


continue_sim = true

dataname = "elong_sphere_11"
datadir = "/home/laigars/sim_data/$dataname"

# elparameters(scale) = Elparameters( # comments are some values that work fine
#  m_use_fraction = false,                     # false
#  m_min_edge_length = 0.7*scale,             # 0.7 * scale
#  m_max_edge_length = 1.5*scale,               # 1.5 * scale
#  m_max_volume_change = 0.1*scale^3,         # 0.1 * scale^3
#  m_min_curvature_multiplier = 1,             # 1
#  m_max_curvature_multiplier = 1,            # 1
#  m_merge_proximity_epsilon = 0.5*scale,     # 0.5 * scale
#  m_proximity_epsilon = 0.00001,             # 0.00001
#  m_perform_improvement = true,              # true
#  m_collision_safety = false,                 # false
#  m_min_triangle_angle = 15,                 # 15
#  m_max_triangle_angle = 120,                # 120
#  m_allow_vertex_movement = false,           # false   ### This is where is a bug
#  m_use_curvature_when_collapsing = false,    # false
#  m_use_curvature_when_splitting = false,    # false
#  m_dt = 1                                   # 1
# )
#par = elparameters(scale




#@load "/home/lai/Dropbox/dokt/code/data/elongation2/data14000.jld2" data
#global points2, faces2, H0 = data[1], data[2], data[3]

H0 = [0., 0., 1.]
mu = 30.

# for mu=30 hist jump @ Bm=3.69
Bm = 0.5
R0 = 1.
#lambda = 10.1
lambda = 1.
gamma = H0[3]^2 * R0 / Bm


reset_vmax = true

t = 0
dt = 0.15
steps = 5000

#points, faces = data[1], data[2]
normals = Normals(points, faces)

max_vs = zeros(3, steps)
mean_vs = zeros(3, steps)

if continue_sim
    reset_vmax = false

    last_file = readdir(datadir)[end]
    global data
    @load "$datadir/$last_file" data
    println("continuing simulation from: $datadir/$last_file")

    global points, faces, t, H0, Bm, v0max = data[1], data[2], data[3], data[4], data[5], data[6]


    # to investigate hysteris region for "elong_sphere_11"
    Bm = 3.8

    global last_step = parse(Int32, last_file[5:9])
    normals = Normals(points, faces)
end


if !isdir("$datadir")
    mkdir("$datadir")
    open("$datadir/_params.txt", "w") do file
        write(file, "H0=$H0\nmu=$mu\nBm=$Bm\neta=$eta\ngamma=$gamma\nlambda=$lambda\nsteps=$steps")
    end
    println("Created new dir: $datadir")
end

for i in 1:steps
    println("------------------------------------------------------------------------------------------------- Step ($i)$(i+last_step)")
    global points, faces, connectivity, normals, all_vs, velocities
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
    mup = mu
    # magnitudes squared of the normal force
    Hn_2 = sum(Hn.^2, dims=1)
    # magnitudes squared of the tangential force
    Ht_2 = sum(Ht.^2, dims=1)

    #tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2
    #tensorn = tensorn * Bm

    # make the force normal to surface // unneeded for Zinch
    #global tensorn = normals .* tensorn

    #zc = SG.Zinchenko2013(points, faces, normals)

    println("velocity:")
    #@time velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)
    #velocities = normals .* velocitiesn_norms' # need to project v_norms on surface

    #@time velocities = make_magvelocities(points, normals, lambda, Bm, mu, Hn_2, Ht_2)

    # dt = 0.1*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
    # if dt > 0.4
    #     dt = 0.4
    # end



    # "_2" has early exit for lambda=1
    @time velocities = make_magvelocities_2(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
    @time velocities = make_Vvecs_conjgrad(normals,faces, points, velocities, 1e-6, 500)

    #passive stabilization
    #velocities = SG.stabilise!(velocities, points, faces, normals, zc)
    #global velocities = make_Vvecs_conjgrad(normals,faces, points, velocities, 1e-6, 120)

    dt = 0.2*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
    if dt < 0.2
        dt = 0.2
    end


    t += dt
    println("dt = $dt")
    println("t = $t")

    points = points + velocities * dt




    # while continue_stab
    #     println("step $i: stabbing .. ")
    #     continue_stab = false
    #
    #     edges = make_edges(faces)
    #     connectivity = make_connectivity(edges)
    #     normals, CDE = make_normals_spline(points, connectivity, edges, normals)
    #
    #     continue_stab = flip_edges(faces, connectivity, points)
    #     points = active_stabilize(points, faces, CDE, connectivity, edges, normals,deltakoef=0.05)
    #
    # end

    #
    do_active = false
    faces, connectivity, do_active = flip_edges(faces, connectivity, points)

    if i % 50 == 0 || do_active
        if do_active
            println("------------------------------------------------------ flipped at step $i")
            edges = make_edges(faces)
            #connectivity = make_connectivity(edges)
            #println("after re-did sum: ", sum(edges))
        end

        #println("OUTSIDE re-did sum: ", sum(edges))
        println("-- doing active / step $i / flipped?: $do_active")
        normals, CDE = make_normals_spline(points, connectivity, edges, normals)
        points = active_stabilize(points, faces, CDE, connectivity, edges, normals,deltakoef=0.05)

    end

    #velocitiesn = velocitiesn'

    #dt = 0.1 * scale / max(sqrt(sum(Vvecs.*Vvecs,2)))
    # ElTopo magic
    #actualdt,points2,faces2 = improvemeshcol(points,faces,points2,par)


    if reset_vmax
        println("Resetting v0max")
        global v0max = maximum(abs.(velocities))
        reset_vmax = false
    end

    vi = maximum(abs.(velocities))
    max_vs[:,i] = [velocities[1, argmax(abs.(velocities[1,:]))],
                    velocities[2, argmax(abs.(velocities[2,:]))],
                    velocities[3, argmax(abs.(velocities[3,:]))]]
    mean_vs[:,i] = [StatsBase.mean(abs.(velocities[1,:])),
                    StatsBase.mean(abs.(velocities[2,:])),
                    StatsBase.mean(abs.(velocities[3,:]))]

    if vi > v0max
        println("updated v0max")
        v0max = vi
    end
    println("Bm = $Bm")
    println("vi = $vi, v0max = $v0max, vi/v0max = $(vi/v0max)")
    println("mean vs: $(mean_vs[:,i])")

    a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])
    println(" --- c/a = $(c/a) , c/b = $(c/b)")

    if i % 1 == 0
        data = [points, faces, t, H0, Bm, v0max]
        println("Finished step $(last_step + i)")
        @save "$datadir/data$(lpad(i + last_step,5,"0")).jld2" data

    end

    if vi/v0max < 0.15
        println("-------------------------------------------------------------------- Increasing Bm at step $i")
        global reset_vmax
        reset_vmax = true
        global Bm
        if Bm < 3.9
            Bm += 0.05
            println("----- new Bm = $Bm")
        end
        # if expand
        #     println("expandasdas")
        #     H0 += [0,0,0.5]
        # else
        #     H0 -= [0,0,0.05]
        # end
        # if H0[3] >= 3.5
        #     global expand = false
        # end
    end


end

data = [max_vs, mean_vs, points, faces]
@save "/home/laigars/sim_data/$(dataname)_v.jld2" data

println("Sim done :)")

scene = Makie.mesh(points', faces', color = :gray, shading = false, visible = true)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 2)


# scene = Makie.mesh(points2_small', faces1',color = :gray, shading = false, visible = true)
# Makie.wireframe!(scene[end][1], color = :blue, linewidth = 2,visible = true)
#
# scene = Makie.mesh(points1_large', faces1',color = :gray, shading = false, visible = true)
# Makie.wireframe!(scene[end][1], color = :red, linewidth = 2,visible = true)


using Plots
using PyPlot
pygui()

fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")

(x, y, z) = [points[i,:] for i in 1:3]
(vx, vy, vz) = [velocities[i,:] for i in 1:3]

ax[:scatter](x,y,z, s=2,color="k")
ax[:quiver](x,y,z,vx,vy,vz, length=30, arrow_length_ratio=0.5)

ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()

#
# using PyPlot
#
# pygui(true)
#
# fig = figure()
# ax = fig[:gca](projection="3d")
#
# N = 10
# x,y,z,u,v,w = [randn(N) for _ in 1:6]
# ax[:quiver](x,y,z, u,v,w)
