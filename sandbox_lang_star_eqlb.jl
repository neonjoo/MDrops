using JLD2
using StatsBase
using LinearAlgebra
using FastGaussQuadrature
using Optim

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")


dir = "star_6"
sourcedir = "/home/laigars/sim_data/$dir"

len = size(readdir(sourcedir),1) - 2

vs = []
ts = zeros(Float64, len)
dSs = []
Vs = zeros(Float64, len)

dira = readdir(sourcedir)

for i in 4:1:len
    #println("$i ")
    f = dira[i]

    #@load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    @load "$sourcedir/$f" data
    points = data[1]
    faces = data[2]
    t = data[3]
    v_phys = data[4]

    normals = Normals(points, faces)
    dS = make_dS(points, faces)
    S = sum(dS)
    global v_phys, dS, points, faces, normals, S

    #println("$i: $(size(points,2)), velo: $(size(v_phys,2))")



    push!(vs, v_phys)
    ts[i] = t
    Vs[i] = make_volume(points, faces, normals)
    push!(dSs, dS)
end

v_avg = zeros(Float64, (3, size(dSs, 1)))

for i = 1:size(dSs,1)
    global vs
    v_phys = vs[i+1]
    dS = dSs[i]
    println(size(v_phys))
    println(size(dS))
    mean_vx = 1/S * sum(v_phys[1,:] .* dS)
    mean_vy = 1/S * sum(v_phys[2,:] .* dS)
    mean_vz = 1/S * sum(v_phys[3,:] .* dS)

    v_avg[:, i] .= mean_vx, mean_vy, mean_vz
end


#
# scene = Makie.mesh(points_img', faces_img', color = :gray, shading = false, visible = true, )
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 1,limits=FRect3D((-2,-2,-4),(4,4,8)))
#
# cam = Makie.cameracontrols(scene)
# cam.upvector[] = (1.0, 0.0, 0)
# update_cam!(scene, cam)
# display(scene)
