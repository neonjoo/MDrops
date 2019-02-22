using LinearAlgebra
using CSV
using JLD2
using Makie
using StatsBase
#using ElTopo


include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry
#include("./SurfaceGeometry/dt20L/src/StabilisationMethods/stabilisationV2.jl")

#include("./SurfaceGeometry/dt20L/src/Properties.jl")

#include("./SurfaceGeometry/dt20L/src/Utils.jl")

include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")


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


H0 = [0, 0, 20]


#H0 = 33 .* [0, 0, 1]
#H0 = [0,0,0]
mu = 7
eta = 1
gamma = 6.9 * 10^-1

continue_sim = false
last_step = 0

w = 2*pi/50
t = 0
dt = 0.01
steps = 20

datadir="/home/laigars/sim_data/"
# typical triangle side length
#scale = 0.1 * 4.9 * 10^-1
#
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
#
#par = elparameters(scale)



global points2 = copy(points)
global faces2 = copy(faces)



if !isdir("$datadir")
    mkdir("$datadir")
    open("$datadir/_params.txt", "w") do file
        write(file, "H0=$H0\nmu=$mu\neta=$eta\ngamma=$gamma\nsteps=$steps")
    end
    println("Created new dir: $datadir")
end

if continue_sim
    last_file = readdir(datadir)[end]
    println(last_file)

    # fix hardcoding
    global data
    @load "$datadir/data02000.jld2" data
    global last_step = 0
end

#@load "/home/lai/Dropbox/dokt/code/data/elongation2/data14000.jld2" data
#global points2, faces2, H0 = data[1], data[2], data[3]

expand = true
reset_vmax = true

normals = Normals(points, faces)


for i in 1:steps
    println("----------------------------------------------------------------------------------------------------- Step $i")
    global points, faces, con, normals
    global t, w, H0
    #global points = copy(points2)
    #global faces = copy(faces2)
    #global expand, reset_vmax

    edges = make_edges(faces)
    con = make_connectivity(edges)
    normals, CDE = make_normals_spline(points, con, edges, normals)

    psi = PotentialSimple(points, faces, mu, H0; normals = normals)
    Ht = HtField(points, faces, psi, normals)
    Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
    Hn = normals .* Hn_norms'

    mup = mu
    # magnitudes squared of the normal force
    Hn_2 = sum(Hn.^2, dims=1)
    # magnitudes squared of the tangential force
    Ht_2 = sum(Ht.^2, dims=1)


    tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2

    # make the force normal to surface
    #tensorn = normals .* tensorn



    zc = SG.Zinchenko2013(points, faces, normals)

    println("velocity:")
    @time velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)

    velocities = normals .* velocitiesn_norms'

    dt = 0.5*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
    println("dt = $dt")

    #global points_ns = points + velocities * dt

    #passive stabilization


    velocitiesn = SG.stabilise!(velocities, points, faces, normals, zc)

    points = points + velocities * dt

    #

    # con0 = copy(con)
    #
    # println("before ", sum(edges))
    # edges0 = copy(edges)
    println("paraboloiding ...")
    #
    #
    do_active = false
    do_active = flip_edges!(faces, con, points)

    # println("after")
    # println(con)
    #
    # println("delta")
    if i % 1 == 0 || do_active
        if do_active
            #println("re-did edges, step $i")
            edges = make_edges(faces)
            #println("after re-did sum: ", sum(edges))
        end

        #println("OUTSIDE re-did sum: ", sum(edges))
        println("doing active / step $i / flipped?: $do_active")
        normals, CDE = make_normals_spline(points, con, edges, normals)
        points = active_stabilize(points, faces, CDE, con, normals,deltakoef=0.1)

    end

    #velocitiesn = make_Vvecs_conjgrad(normals,faces, points, velocitiesn, 1e-6, 120);

    #velocitiesn = velocitiesn'

    #dt = 0.1 * scale / max(sqrt(sum(Vvecs.*Vvecs,2)))



    # ElTopo magic
    #actualdt,points2,faces2 = improvemeshcol(points,faces,points2,par)
    #
    # if i % 1 == 0
    #     data = [points2, faces2, (H0)]
    #     println("Finished step $(last_step + i)")
    #     @save "$datadir/data$(lpad(i + last_step,5,"0")).jld2" data
    #
    # end
    #
    #
    # if reset_vmax
    #     println("Resetting v0max")
    #     global v0max = maximum(abs.(velocitiesn))
    #     reset_vmax = false
    # end
    #
    # vi = maximum(abs.(velocitiesn))
    #
    # println("vi = $vi, v0max = $v0max, vi/v0max = $(vi/v0max)")
    #
    # if vi/v0max < 0.05
    #     println("Reached eqlb with vi=$vi")
    #     global expand, reset_vmax
    #     reset_vmax = true
    #
    #     if expand
    #         println("expandasdas")
    #         H0 += [0,0,0.5]
    #     else
    #         H0 -= [0,0,0.05]
    #     end
    #     if H0[3] >= 3.5
    #         global expand = false
    #     end
    #end


end



#
# Hn_teor = zeros(size(points,2))
# Ht_teor = zeros(3, size(points,2))
#
# for i in 1:size(points,2)
#     r = norm(points[:,i])
#     r_xy = norm(points[1:2,i])
#     (x,y,z) = [points[j,i] for j in 1:3]
#     Hn_teor[i] = 3 * norm(H0) * points[3,i] / r / (mu+2)
#     Ht_teor[1,i] = - 3 * norm(H0) * r_xy / r * x*z /r/r_xy/(mu+2)
#     Ht_teor[2,i] = - 3 * norm(H0) * r_xy / r * y*z / r/r_xy/(mu+2)
#     Ht_teor[3,i] = 3 * norm(H0) * r_xy / r * r_xy / r / (mu+2)
# end



scene = Makie.mesh(points', faces', color = :gray, shading = false,visible = false)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 1)
scene = Makie.mesh!(points_ns', faces',color = :gray, shading = false,visible = false)
Makie.wireframe!(scene[end][1], color = :blue, linewidth = 1,visible = true)
