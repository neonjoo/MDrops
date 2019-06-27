using Makie
using JLD2
using FileIO


dir = "elong_sphere_zinch4"
sourcedir = "/home/laigars/sim_data/$dir"
outdir="/home/laigars/sim_data/pics/$dir"
len = size(readdir(sourcedir),1) - 1

if !isdir("$outdir")
    #mkdir("$outdir")
    mkpath("$outdir")
end

ratios = []
steps = []
for i in 5:100:50000
    #i =
    @load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    println("step $i")
    points_img = data[1]
    faces_img = data[2]

    scene = Makie.mesh(points_img', faces_img', color = :gray, shading = false, visible = true)
    Makie.wireframe!(scene[end][1], color = :black, linewidth = 1,limits=FRect3D((-1,-1,-3),(3,3,6)))



    # global ratios, steps
    # ratio = (maximum(points2[3,:]) - minimum(points2[3,:])) / (maximum(points2[1,:]) - minimum(points2[1,:]))
    # push!(ratios, ratio)
    # push!(steps, i)

    # text!(
    # "$i",
    # position = (-2, 0,5),
    # textsize = 1)
    save("$outdir/$(lpad(i,5,"0")).png", scene)

end


#
# scene = Makie.mesh(points_img', faces_img', color = :gray, shading = false, visible = true, )
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 1,limits=FRect3D((-2,-2,-4),(4,4,8)))
#
# cam = Makie.cameracontrols(scene)
# cam.upvector[] = (1.0, 0.0, 0)
# update_cam!(scene, cam)
# display(scene)
