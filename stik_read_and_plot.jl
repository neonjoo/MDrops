using Pkg
pkg"activate ."
pkg"resolve"

#using Makie
using FileIO
using JLD2
using StatsBase
using Optim

datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_manymoreN_adaptive_dt/"
#datadir2="/home/andris/sim_data/2019-03-15/1

#file_no = 1180
#last_file = readdir(datadir)[file_no]
last_file = readdir(datadir)[end-3]
#last_file = "data00114.jld2"
println()
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

#scene = Makie.mesh(points', faces',color = :lightgray, shading = true, visible = true)
#Makie.wireframe!(scene[end][1], color = :black, linewidth = 0.7,visible = true)#, limits=FRect3D((-5,-5,-5),(10,10,10)))
#text!("$(last_file[5:end-5])", position = (0, 0, 4), textsize = 0.6, rotation=2.)


#display(scene)
#cam = Makie.cameracontrols(scene)
#cam.upvector[] = (-1., -1., 1.)
#update_cam!(scene, cam)
# scene.center = false
# scene
# FileIO.save("plot.png", scene)
#readline(stdin)


#using Plots
ENV["MPLBACKEND"]="tkagg"
using PyPlot
pygui(true)

fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")

(x, y, z) = [points[i,:] for i in 1:3]

ax[:scatter](x,y,z, s=2,color="k")

ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()
