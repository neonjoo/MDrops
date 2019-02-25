using Makie
using FileIO
using JLD2

datadir="/home/andris/mydatadirst_weekend/"

files = readdir(datadir)
for file in files
    @load "$datadir/$file" data

    points2, faces2, dt = data[1], data[2], data[3]
    faces2 = Array{Int64,2}(faces2)


    edges2 = make_edges(faces2)
    connectivity = make_connectivity(edges2)

    minvalence = minimum(sum(x->x!=0,connectivity,dims=1))

    println("file = $file"," minvalence = $minvalence")

end






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
