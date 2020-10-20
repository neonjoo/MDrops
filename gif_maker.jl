using Makie
using JLD2
using FileIO
using StatsBase

dir = "star_5"
sourcedir = "/home/laigars/sim_data/$dir"
outdir="/home/laigars/sim_data/pics/$dir"

#sourcedir = "/mnt/big_data/shared_folder/plotation/$dir"
#outdir="/mnt/big_data/shared_folder/plotation/$dir"
len = size(readdir(sourcedir),1) - 1

if !isdir("$outdir")
    #mkdir("$outdir")
    mkpath("$outdir")
end

ratios = []
steps = []
for f in readdir(sourcedir)[1747:1:end]
    println("f = $f")
    #@load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    @load "$sourcedir/$f" data
    points = data[1]
    faces = data[2]

    mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                    StatsBase.mean(points[2,:]),
                    StatsBase.mean(points[3,:]))

    points = points .- [mean_x, mean_y, mean_z]
    scene = Makie.mesh(points', faces', color = :gray, shading = false, visible = true)
    Makie.wireframe!(scene[end][1], color = :black, linewidth = 1,limits=FRect3D((-1,-1,-2),(2,2,4)))


    # global ratios, steps
    # ratio = (maximum(points2[3,:]) - minimum(points2[3,:])) / (maximum(points2[1,:]) - minimum(points2[1,:]))
    # push!(ratios, ratio)
    # push!(steps, i)

    # text!(
    # "$i",
    # position = (-2, 0,5),
    # textsize = 1)

    #display(scene)
    save("$outdir/$f.png", scene)

end


#
# scene = Makie.mesh(points_img', faces_img', color = :gray, shading = false, visible = true, )
# Makie.wireframe!(scene[end][1], color = :black, linewidth = 1,limits=FRect3D((-2,-2,-4),(4,4,8)))
#
# cam = Makie.cameracontrols(scene)
# cam.upvector[] = (1.0, 0.0, 0)
# update_cam!(scene, cam)
# display(scene)
