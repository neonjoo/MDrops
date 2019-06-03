using LinearAlgebra
using CSV
using JLD2
using Makie
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

#points_csv= CSV.read("/home/lai/Dropbox/dokt/code/matlab/vertices.csv", header=0)
#faces_csv = CSV.read("/home/lai/Dropbox/dokt/code/matlab/triangles.csv", header=0)


points_csv= CSV.read("./meshes/points_sphere.csv", header=0)
faces_csv = CSV.read("./meshes/faces_sphere.csv", header=0)

println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')


points0 = copy(points)
faces0 = copy(faces)

continue_sim = true
last_step = 0

dataname = "elong_sphere_zinch3"
datadir = "/home/joo/sim_data/$dataname"
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
#par = elparameters(scale)




#@load "/home/lai/Dropbox/dokt/code/data/elongation2/data14000.jld2" data
#global points2, faces2, H0 = data[1], data[2], data[3]

H0 = [0., 0., 1.]
mu = 3.
lambda = 1.
Bm = 3.
R0 = 1.
eta = 1.
gamma = H0[3]^2 * R0 / Bm


expand = true
reset_vmax = false

t = 0
dt = 0.04
steps = 500

#points, faces = data[1], data[2]
#normals = Normals(points, faces)

max_vs = zeros(3, steps)
mean_vs = zeros(3, steps)

if continue_sim
    last_file = readdir(datadir)[end]

    # fix hardcoding
    # fix v_max overwrite upon continuing

    global data
    @load "$datadir/$last_file" data
    println("continuing simulation from: $datadir/$last_file")

    global points, faces, t, H0, Bm, v0max = data[1], data[2], data[3], data[4], data[5], data[6]
    global last_step = 1000
end


if !isdir("$datadir")
    mkdir("$datadir")
    open("$datadir/_params.txt", "w") do file
        write(file, "H0=$H0\nmu=$mu\neta=$eta\ngamma=$gamma\nsteps=$steps")
    end
    println("Created new dir: $datadir")
end

for i in 1:steps
    println("----------------------------------------------------------------------------------------------------- Step $i")
    global points, faces, connectivity, normals, all_vs
    global t, H0
    #global points = copy(points2)
    #global faces = copy(faces2)
    #global expand, reset_vmax

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
    global Hn_2 = sum(Hn.^2, dims=1)
    # magnitudes squared of the tangential force
    global Ht_2 = sum(Ht.^2, dims=1)


    tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2
    tensorn = tensorn * Bm

    # make the force normal to surface // unneeded for Zinch
    #global tensorn = normals .* tensorn

    #zc = SG.Zinchenko2013(points, faces, normals)

    println("velocity:")
    @time velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)
    global velocities = normals .* velocitiesn_norms' # need to project v_norms on surface

    #@time velocities = make_magvelocities(points, normals, lambda, Bm, mu, Hn_2, Ht_2)

    # dt = 0.1*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
    # if dt > 0.4
    #     dt = 0.4
    # end
    println("dt = $dt")
    println("t = $t")
    t += dt

    #passive stabilization

    #velocities = SG.stabilise!(velocities, points, faces, normals, zc)
    #global velocities = make_Vvecs_conjgrad(normals,faces, points, velocities, 1e-6, 120)
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

    if i % 25 == 0 || do_active
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
    mean_vs[:,i] = [mean(abs.(velocities[1,:])),
                    mean(abs.(velocities[2,:])),
                    mean(abs.(velocities[3,:]))]

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

    if vi/v0max < 0.10
        println("-------------------------------------------------------------------- Increasing Bm at step $i")
        global reset_vmax
        reset_vmax = true
        global Bm
        Bm += 2
        println("----- new Bm = $Bm")
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
@save "/home/joo/sim_data/$(dataname)_v4.jld2" data

println("Sim done :)")

# scene = Makie.mesh(points', faces', color = :gray, shading = false, visible = true)
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 2)

#
# scene = Makie.mesh(points2_small', faces1',color = :gray, shading = false, visible = true)
# Makie.wireframe!(scene[end][1], color = :blue, linewidth = 2,visible = true)
#
# scene = Makie.mesh(points1_large', faces1',color = :gray, shading = false, visible = true)
# Makie.wireframe!(scene[end][1], color = :red, linewidth = 2,visible = true)


using Plots
using PyPlot

pyplot()
pygui()

fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")
(x, y, z) = [points[i,:] for i in 1:3]
(vnx, vny, vnz) = [velocities[i,:] for i in 1:3]
ax[:scatter](x,y,z, s=2,color="k")
ax[:quiver](x,y,z,vnx,vny,vnz, length=3, arrow_length_ratio=0.5)
ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()
