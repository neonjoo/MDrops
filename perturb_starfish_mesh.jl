#cd("/home/andris/MDrops/")

using Pkg

pkg"activate ."
pkg"resolve"

using JLD2
using StatsBase
using LinearAlgebra
using FastGaussQuadrature
using Optim
using Distributed
using Dates

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")

cs = [0.67381703, 0.65155044, 0.62959378, 0.60794706, 0.58661027, 0.56558342, 0.5448665 , 0.52445952, 0.50436247]
bms = [17, 18, 19, 20, 21, 22, 23, 24, 25]

for (idx, bm) in enumerate(bms)
	global points, faces
	println("##########################################################################################, Bm=$bm")


	@load "./meshes/two_axis_ellipsoid_Bm_$bm.jld2" data

	points, faces = data[1], data[2]
	points = Array{Float64}(points)
	faces = Array{Int64}(faces)

	# the volume change is 0.05%

	k = 0.03
	n = 2
	ts = atan.(points[3,:], points[1,:])
	rs = sqrt.(points[1,:].^2 + points[3,:].^2)
	perturb = k .* cos.(n*ts)

	global pertpoints = copy(points)
	pertpoints[3,:] = pertpoints[3,:] .+ perturb .* sin.(ts) .* rs
	pertpoints[1,:] = pertpoints[1,:] .+ perturb .* cos.(ts) .* rs

#
#scene = Makie.mesh(pertpoints', faces',color = :lightgray, shading = false, visible = true)
#Makie.wireframe!(scene[end][1], color = :black, linewidth = 0.7,visible = true)#, limits=FRect3D((-5,-5,-5),(10,10,10)))
#
# Makie.mesh!(points', faces',color = :green, shading = false, visible = false)
# Makie.wireframe!(scene[end][1], color = :red, linewidth = 0.7,visible = true)#, limits=FR

	data = [pertpoints, faces]
	#@save "./meshes/perturbed/ellipsoid_n_$(n)_Bm_$(bm).jld2" data

end
