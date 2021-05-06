using Pkg
pkg"activate ."
pkg"resolve"

using LinearAlgebra
#using CSV
using JLD2
#using Makie
using StatsBase
using Optim
using FastGaussQuadrature
using Distributed

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./stabilization.jl")
#include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
#include("./sandbox_lang.jl")

#points_csv= CSV.read("./meshes/points_critical_hyst_2_21.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_critical_hyst_2_21.csv", header=0)
#fields = CSV.read("/home/laigars/sim_data/field.csv", header=0)[1] * 10 # mT -> Oe
#times = CSV.read("/home/laigars/sim_data/time.csv", header=0)[1]
#points_csv= CSV.read("./meshes/points_ellipse_fewN.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_ellipse_fewN.csv", header=0)

# points = convert(Array, points_csv)
# faces = convert(Array, faces_csv)
points, faces = expand_icosamesh(R=1, depth=3)

#@load "./meshes/faces_critical_hyst_2_21.jld2" faces
points = Array{Float64}(points)
faces = Array{Int64}(faces)
#points = points

println("Running on $(Threads.nthreads()) threads")

continue_sim = false

dataname = "test_spherical"
datadir = "/home/laigars/sim_data/$dataname"

H0 = [0., 0., 1.]
mu = 10

# Bm_crit = 3.68423 pie mu=30
Bm = 30 ################################################ zemāk iespējams loado citu
#R0 = 21.5 * 100/480 * 1e-4 # um to cm for cgs
R0 = 1
lambda = 7.6
gamma = H0[3]^2 * R0 / Bm
#gamma = 8.2 * 1e-4
#gamma = 7.7 * 1e-4 # from fitted exp data with mu=34
w = 0

reset_vmax = true
last_step = 0
t = 0
time_k = 0.05012312572132173
dt = 2pi / w * time_k
dt = 0.01
steps = 1
epsilon = 0.05

normals = Normals(points, faces)

max_vs = zeros(3, steps)
mean_vs = zeros(3, steps)
max_abs_v = zeros(1, steps)

if continue_sim
    reset_vmax = false

    last_file = readdir(datadir)[end]
    global data
    println("continuing simulation from: $datadir/$last_file")
    @load "$datadir/$last_file" data
    # data = [points, faces, t, velocities_phys, H0, Bm]
    global points, faces, t, velocities_phys, H0, Bm = data[1], data[2], data[3], data[4], data[5], data[6]
    global last_step = parse(Int32, last_file[5:9])
    println("last step: $last_step")
    normals = Normals(points, faces)
end


println("Loaded mesh; nodes = $(size(points,2))")

if !isdir("$datadir")
    mkdir("$datadir")
    println("Created new dir: $datadir")
    open("$datadir/aa_params.txt", "w") do file
        write(file, "H0=$H0\nmu=$mu\nBm=$Bm\nlambda=$lambda\nsteps=$steps\ndt=$dt\nw=$w\n")
    end
    cp("main_lan.jl", "$datadir/aa_source_code.jl")
end

previous_i_when_flip = 0

for i in 1:steps
    println("---------------------------------Number of nodes: $(size(points,2))----------------------- Step ($i)$(i+last_step)")
    global points, faces, connectivity, normals, all_vs, velocities, velocities_phys
    global t, H0, epsilon
    global max_abs_v, max_v_avg
    global previous_i_when_flip
    edges = make_edges(faces)
    connectivity = make_connectivity(edges)
    normals, CDE, AB = make_normals_parab(points, connectivity, normals)
    neighbor_faces = make_neighbor_faces(faces)

    psi = PotentialSimple_par(points, faces, normals, mu, H0)
    Ht = HtField_par(points, faces, psi, normals)
    print("Normal field calculation: ")
    @time Hn_norms = NormalFieldCurrent_par(points, faces, normals, Ht, mu, H0)
    Hn = normals .* Hn_norms'

    Hn_2 = sum(Hn.^2, dims=1)
    Ht_2 = sum(Ht.^2, dims=1)

    #println("Bm = $Bm")

    @time velocities_phys = make_magvelocities_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
    @time velocities = make_Vvecs_conjgrad(normals,faces, points, velocities_phys, 1e-6, 500)

    #dt = 0.05*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
     #if dt < 0.2
    #     dt = 0.2
    # end
    #println("max v: $(maximum(abs.(velocities))),   min v: $(minimum(abs.(velocities)))")
    t += dt
    println("---- t = $t, dt = $dt ------")

    points = points + velocities * dt
    normals, CDE, AB = make_normals_parab(points, connectivity, normals)

    H0 = [sin(w*t), 0., cos(w*t)]

    cutoff_crit = 0.4 # 0.5 for sqrt(dS), 0.55 for max triangle edge length
    minN_triangles_to_split = 5

    marked_faces  = mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = cutoff_crit)
    if sum(marked_faces) > minN_triangles_to_split
        println("-----------------------------------")
        println("----------Adding mesh points-------")
        println("    V-E+F = ", size(points,2)-size(edges,2)+size(faces,2))
        println("    number of points: ", size(points,2))
        println("    number of faces: ", size(faces,2))
        println("    number of edges: ", size(edges,2))
        println("-----------------------------------")

        points_new, faces_new = add_points(points, faces,normals, edges, CDE; cutoff_crit = cutoff_crit)
        edges_new = make_edges(faces_new)
        connectivity_new = make_connectivity(edges_new)

        println("New V-E+F = ", size(points_new,2)-size(edges_new,2)+size(faces_new,2))
        println("New number of points: ", size(points_new,2))
        println("New number of faces: ", size(faces_new,2))
        println("New number of edges: ", size(edges_new,2))
        println("-----------------------------------")
        println("active stabbing after adding points")
        println("------flipping edges first---------")
        faces_new, connectivity_new, do_active = flip_edges(faces_new, connectivity_new, points_new)
        edges_new = make_edges(faces_new)
        println("-- flipped?: $do_active")
        println("---- active stabbing first --------")
        points_new = active_stabilize_old_surface(points,CDE,normals,points_new, faces_new, connectivity_new, edges_new)
        println("------flipping edges second---------")
        faces_new, connectivity_new, do_active = flip_edges(faces_new, connectivity_new, points_new)
        edges_new = make_edges(faces_new)
        println("-- flipped?: $do_active")
        println("---- active stabbing second --------")
        points_new = active_stabilize_old_surface(points,CDE,normals,points_new, faces_new, connectivity_new, edges_new)

        points, faces, edges, connectivity = points_new, faces_new, edges_new, connectivity_new
        normals = Normals(points, faces)
        normals, CDE, AB = make_normals_parab(points, connectivity, normals)
        println("New normals pointing out? ", all(sum(normals .* points,dims=1).>0))
        println("-----------------------------------")
        println("---------- Points added -----------")
        println("-----------------------------------")

    else # stabilize regularly if havent added new faces
        do_active = false

        faces, connectivity, do_active = flip_edges(faces, connectivity, points)

        if do_active
            if i - previous_i_when_flip > 5
                println("doing active because flipped")
                edges = make_edges(faces)
                points = active_stabilize(points, faces, CDE, connectivity, edges, normals)
            else
                println("flipped; not doing active")
            end
            previous_i_when_flip = i
        end

        if i % 100 == 0 && i > 2
            println("doing active every 100th time step")
            points = active_stabilize(points, faces, CDE, connectivity, edges, normals)
        end
    end

    #dt = 0.1 * scale / max(sqrt(sum(Vvecs.*Vvecs,2)))
    # ElTopo magic
    #actualdt,points2,faces2 = improvemeshcol(points,faces,points2,par)


    #println("Bm = $Bm")
    #println("vi = $vi, max_v_avg = $max_v_avg, v0max = $v0max, max_v_avg/v0max = $(max_v_avg/v0max)")
    #println("mean vs: $(mean_vs[:,i])")

    a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])
    println(" --- c/a = $(c/a) , c/b = $(c/b)")

    if i % 1 == 0
        data = [points, faces, t, velocities_phys, H0, Bm]
        @save "$datadir/data$(lpad(i + last_step,5,"0")).jld2" data
    end

end # end simulation iterations

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


function v_x(points)
    x, y, z = points[1,:], points[2,:], points[3,:]
    mu = 10
    lambda = 7.6
    Bm = 30

    return 3*Bm*(mu-1)^2 .* x .* (1 .+ 3*z.^2) / ( 5 * lambda * (2+mu)^2)
end

function v_z(points)
    x, y, z = points[1,:], points[2,:], points[3,:]
    mu = 10
    lambda = 7.6
    Bm = 30

    return 3*Bm*(mu-1)^2 .* z .* (1 .- 6*(x.^2 .+ y.^2) .- 3*z.^2 ) / ( 5 * lambda * (2+mu)^2)
end
