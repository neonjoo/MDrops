using Makie
#using FileIO
using JLD2


datadir="/home/laigars/sim_data/random_long3"
#datadir2="/home/laigars/sim_data/random_long3"

last_file = readdir(datadir)[end]
last_file2 = readdir(datadir2)[end]

last_file = "data00084.jld2"
last_file2 = "data00045.jld2"

println("first: $last_file")
println("second: $last_file2")


@load "$datadir/$last_file" data
#@load "$datadir/data0005.jld2" data

points2, faces2 = data[1], data[2]
faces2 = Array{Int64,2}(faces2)
println("plotting")
println(size(points2))
println(size(faces2))
scene = Makie.mesh(points2', faces2',color = :grey, shading = false,visible = true)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 1)


@load "$datadir2/$last_file2" data

#points2, faces2 = data[1], data[2]
#faces2 = Array{Int64,2}(faces2)
#scene = Makie.mesh!(points2', faces2',color = :white, shading = false,visible = false)
#Makie.wireframe!(scene[end][1], color = :blue, linewidth = 1)




display(scene)
# cam = Makie.cameracontrols(scene)
# cam.upvector[] = (0.0, 0.0, 1.0)
# update_cam!(scene, cam)
# scene.center = false
# scene
# FileIO.save("plot.png", scene)
readline(stdin)
