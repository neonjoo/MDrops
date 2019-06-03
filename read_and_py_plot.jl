using JLD2

datadir="/home/joo/sim_data/elong_sphere_zinch3/"
#datadir2="/home/andris/sim_data/2019-03-15/1"

last_file = readdir(datadir)[end]
println(last_file)

# fix hardcoding
global data
@load "$datadir/$last_file" data
#@load "$datadir/data00031.jld2" data


points, faces = data[1], data[2]


println("plotting")

using Plots
using PyPlot

pyplot()
pygui()

fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")
(x, y, z) = [points[i,:] for i in 1:3]
#(vnx, vny, vnz) = [velocities[i,:] for i in 1:3]
ax[:scatter](x,y,z, s=2,color="k")
#ax[:quiver](x,y,z,vnx,vny,vnz, length=20, arrow_length_ratio=0.5)
ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()

println("plotted")
readline(stdin)
