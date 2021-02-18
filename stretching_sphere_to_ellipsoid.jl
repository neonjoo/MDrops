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
	global points, faces, connectivity, edges
	cc = cs[idx]
	println("########################################################################################## c = $cc, Bm=$bm")
	points, faces = expand_icosamesh(R=1, depth=3)
	points = Array{Float64}(points)
	faces = Array{Int64}(faces)

	edges = make_edges(faces)
	neighbor_faces = make_neighbor_faces(faces)
	connectivity = make_connectivity(edges)
	normals = Normals(points, faces)
	normals, CDE, AB = make_normals_parab(points, connectivity, normals; eps = 10^-8)

	cutoff_crit = 0.2
	a = 1/sqrt(3)#0.6^(1/3)
	b = 1/sqrt(3)#0.6^(1/3)
	c = 1/(a*b)

	b = cc
	a = 1/sqrt(b)
	c = 1/sqrt(b)

	iters = 10
	aiters = range(1,a,length=iters)
	biters = range(1,b,length=iters)
	citers = range(1,c,length=iters)

	for i in 2:iters
		println("--------------------------------------------------------------------------------- $i")
		global points, faces, connectivity, edges
		points[1,:] .*= aiters[i]/aiters[i-1]
		points[2,:] .*= biters[i]/biters[i-1]
		points[3,:] .*= citers[i]/citers[i-1]

		# stabilize
		neighbor_faces = make_neighbor_faces(faces)
		normals = Normals(points, faces)
		println("New first approx normals pointing out? ", all(sum(normals .* points,dims=1).>0))
		normals, CDE, AB = make_normals_parab(points, connectivity, normals; eps = 10^-8)

		marked_faces  = mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = cutoff_crit)

		println(sum(marked_faces))
		if sum(marked_faces) >= 1
			global points_new, faces_new = add_points(points, faces,normals, edges, CDE; cutoff_crit = cutoff_crit)
	        global edges_new = make_edges(faces_new)
	        global connectivity_new = make_connectivity(edges_new)
			global marked_points = make_marked_points(points, faces, points_new, marked_faces)

			global do_active = true
			global oldE = Inf #10.^14
			println("locally relaxing added points")
			while do_active
			    points_new, maxgrad = relax_after_split_weights(marked_points, points_new, faces_new, connectivity_new, points, normals, CDE,cutoff_crit; trianw = 5, trianpow = 2, edgepow = 2)

			    points_old = copy(points_new)

			    E = make_simple_E_weights(points_new,faces_new,points, normals, CDE,cutoff_crit; trianw = 5, trianpow = 2, edgepow = 2)
			    println("Mesh energy ",E)
			    if E >= oldE
			        faces_new = faces_old
			        points_new = points_old
			        break
			    end

			    global faces_old = copy(faces_new) # dunno why this has to be global
			    global faces_new, connectivity_new, do_active = flip_edges(faces_new, connectivity_new, points_new)


			    oldE = E
			end

		    points, faces, edges, connectivity = points_new, faces_new, edges_new, connectivity_new
	        normals = Normals(points, faces)
	        println("New first approx normals pointing out? ", all(sum(normals .* points,dims=1).>0))
	        normals, CDE, AB = make_normals_parab(points, connectivity, normals; eps = 10^-8)
	        #normals, CDE = make_normals_spline(points, connectivity, edges, normals)
	        println("New normals pointing out? ", all(sum(normals .* points,dims=1).>0))
	        println("-----------------------------------")
	        println("---------- Points added -----------")
	        println("-----------------------------------")
		end
	end

	data = [points, faces]
	@save "./meshes/two_axis_ellipsoid_Bm_$bm.jld2" data
#%%

	# k = 0.03
	# n = 5
	# ts = atan.(points[3,:], points[1,:])
	# rs = sqrt.(points[1,:].^2 + points[3,:].^2)
	# #ts = range(0, stop=2pi, step=0.01)
	# perturb = k .* cos.(n*ts)
	# #perturb = abs.(perturb)
	#
	# pertpoints = copy(points)
	# pertpoints[3,:] = pertpoints[3,:] .+ perturb .* sin.(ts) .* rs
	# pertpoints[1,:] = pertpoints[1,:] .+ perturb .* cos.(ts) .* rs
	#
	# #scene = Makie.mesh(pertpoints', faces',color = :lightgray, shading = false, visible = true)
	# #Makie.wireframe!(scene[end][1], color = :black, linewidth = 0.7,visible = true)#, limits=FRect3D((-5,-5,-5),(10,10,10)))
	#
	# data = [pertpoints, faces]
	# @save "./meshes/perturbed/ellipsoid_n_$($n)_Bm_$(bm).jld2" data

end
