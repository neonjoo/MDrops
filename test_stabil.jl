using Makie
using FileIO
using JLD2
using Optim
using StatsBase
include("./mesh_functions.jl")


ppoints, pfaces = expand_icosamesh(;R=1,depth=7)
ppoints = Array{Float64}(ppoints)
ppoints[3,:] *= 2.5 # make an ellipse
ppoints0 = copy(ppoints)
pfaces = Array{Int64}(faces)

datadir="/home/andris/mydatadirst_weekend/"

file = readdir(datadir)[48]
@load "$datadir/$file" data
points, faces, dt = data[1], data[2], data[3]
faces = Array{Int64,2}(faces)
edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)

for i = 1:size(points,2)
    global ppoints
    ppoints[:,i] = project_on_drop(points,CDE,normals,ppoints[:,i])
end


scene = Makie.mesh(ppoints', pfaces',color = :white, shading = false,visible = false)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 1)
scene = Makie.mesh!(points', faces',color = :gray, shading = false,visible = true)
Makie.wireframe!(scene[end][1], color = :blue, linewidth = 1,visible = true)



#
# files = readdir(datadir)
# for file in files
#     @load "$datadir/$file" data
#
#     points2, faces2, dt = data[1], data[2], data[3]
#     faces2 = Array{Int64,2}(faces2)
#
#
#     edges2 = make_edges(faces2)
#     connectivity = make_connectivity(edges2)
#
#     minvalence = minimum(sum(x->x!=0,connectivity,dims=1))
#
#     println("file = $file"," minvalence = $minvalence")
#
# end
# 48 ir pēdējais ar minvalenci 5






#
# # fix hardcoding
# global data
# @load "$datadir/$last_file" data
# #@load "$datadir/data00050.jld2" data
#
# points2, faces2 = data[1], data[2]
# faces2 = Array{Int64,2}(faces2)
#
# println("plotting")
# println(size(points2))
# println(size(faces2))
# scene = Makie.mesh(points2', faces2',color = :white, shading = false,visible = true)
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 1)
# display(scene)
# cam = Makie.cameracontrols(scene)
# cam.upvector[] = (0.0, 0.0, 1.0)
# update_cam!(scene, cam)
# scene.center = false
# scene
# FileIO.save("plot.png", scene)
# readline(stdin)
