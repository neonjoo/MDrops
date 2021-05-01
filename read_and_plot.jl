using Pkg
pkg"activate ."
pkg"resolve"

using Makie
using FileIO
using JLD2
using StatsBase
using Optim


#datadir = "/home/laigars/sim_data/star_8/"
datadir = "/home/laigars/mdrops/meshes/perturbed/starfish_2/"
#datadir = "/mnt/hpc/sim_data/perturbed_real_6_15/"
#datadir="/mnt/hpc/sim_data/rotation_lambd10.0_mu30.0_w1.5_Bm50.0"
datadir="/mnt/hpc/sim_data/perturbed_2_21_fastfield/"


last_file = readdir(datadir)[29] # end-3 if there are aux files (source, speeds);  end if there aren't
#last_file = "ellipsoid_n_5_Bm_17.jld2"

println()
println(datadir)
println(last_file)
# fix hardcoding
global data
@load "$datadir/$last_file" data
#@load "$datadir/data00031.jld2" data


points, faces = data[1], data[2]
faces = Array{Int64,2}(faces)

mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                StatsBase.mean(points[2,:]),
                StatsBase.mean(points[3,:]))

points = points .- [mean_x, mean_y, mean_z]


function f(ab::Array{Float64,1})
    return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
            points[3,:].^2/ab[2]^2 .- 1).^2)
end

x0 = [0.99, 1.01]
res = Optim.optimize(f,x0)
b = Optim.minimizer(res)[1]
a = Optim.minimizer(res)[2]

println("a = $a, b = $b")
# @load "$datadir2/$last_file" data
# points3, faces3 = data[1], data[2]
# faces3 = Array{Int64,2}(faces3)
#
# println("plotting")
# println(size(points2))
# println(size(faces2))

n = 20

θ = [0;(0.5:n-0.5)/n;1]
φ = [(0:2n-2)*2/(2n-1);2]
x = [b*cospi(φ)*sinpi(θ) for θ in θ, φ in φ]
y = [b*sinpi(φ)*sinpi(θ) for θ in θ, φ in φ]
z = [a*cospi(θ) for θ in θ, φ in φ]


#Makie.surface(x, y, z, color=:black, shading=false, transparency=true)

scene = Makie.mesh(points', faces',color = :lightgray, shading = false, visible = true)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 0.7,visible = true)#, limits=FRect3D((-5,-5,-5),(10,10,10)))
#text!("$(last_file[5:end-5])", position = (0, 0, 4), textsize = 0.6, rotation=2.)


#display(scene)
#cam = Makie.cameracontrols(scene)
#cam.upvector[] = (-1., -1., 1.)
#update_cam!(scene, cam)
# scene.center = false
# scene
# FileIO.save("plot.png", scene)
#readline(stdin)
