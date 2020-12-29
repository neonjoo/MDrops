cd("/home/andris/MDrops/")

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
using Roots
using Random
rng = MersenneTwister(1234)


include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")

points, faces = expand_icosamesh(R=1, depth=3)

#@load "./meshes/faces_critical_hyst_2_21.jld2" faces
points = Array{Float64}(points)
faces = Array{Int64}(faces)

println("Loaded mesh; nodes = $(size(points,2))")

continue_sim = true

dataname = "fast_field_test14_newmesh"
datadir = "/home/andris/sim_data/$dataname"
println("Simulation: $dataname")

H0 = [0., 0., 1.]
w = Inf
mu = 10.
Bm = 25.
lambda = 1.
t = 0.


dt = 0.05
steps = 3500000
last_step = 0
cutoff_crit = 0.2 # for triangle size

if continue_sim && isdir("$datadir")
    last_file = readdir(datadir)[end]
    global data
    println("continuing simulation from: $datadir/$last_file")
    @load "$datadir/$last_file" data

    global points, faces, t = data[1], data[2], data[3]
    global lambda, mu, w, Bm = data[end-3], data[end-2], data[end-1], data[end]
    if w != Inf
        global H0 = [sin(w*t), 0., cos(w*t)]
    end

    global last_step = parse(Int32, last_file[5:9])
    println("last step: $last_step")

	cp("main_hpc_stik.jl", "$datadir/aa_source_code.jl"; force=true)
end

if !isdir("$datadir")
    mkdir("$datadir")
    println("Created new dir: $datadir")
    open("$datadir/aa_params.txt", "w") do file
        write(file, "H0=$H0\nmu=$mu\nBm=$Bm\nlambda=$lambda\nsteps=$steps\ndt=$dt\nw=$w\n")
    end

	cp("main_stik.jl", "$datadir/aa_source_code.jl"; force=true)
end




previous_i_when_flip = -1000
previous_i_when_split = -1000
println("Running on $(Threads.nthreads()) threads")
steps = steps - last_step
for i in 1:steps
    if t > 150.
        #println("passed t = 15")
    	break
    end
    println("----------------------------------------------------------------------")
    println("----- Number of points: $(size(points,2)) ---------- Step ($i)$(i+last_step)--- t = $(t)-------")
	println("----------$(Dates.format(now(), "yyyy-mm-dd;  HH:MM:SS"))------------")
    println("----------------------------------------------------------------------")
    global points, faces, connectivity, normals, velocities, velocities_n, neighbor_faces, edges, CDE
    global dt, t, H0, V0
    global previous_i_when_flip, previous_i_when_split, cutoff_crit
	global rng
    edges = make_edges(faces)
    neighbor_faces = make_neighbor_faces(faces)
    connectivity = make_connectivity(edges)
	if i == 1
		normals = Normals(points, faces)
	end
    normals, CDE, AB = make_normals_parab(points, connectivity, normals; eps = 10^-8)
	if i == 1
		V0 = make_volume(points,faces, normals)
		println("start volume = ",V0)
	end

	if w != Inf
	    psi = PotentialSimple_par(points, faces, normals, mu, H0)
	    Ht = HtField_par(points, faces, psi, normals)
	    Hn_norms = NormalFieldCurrent_par(points, faces, normals, Ht, mu, H0)
	    Hn = normals .* Hn_norms'

	    # magnitudes squared of the normal force
	    Hn_2 = sum(Hn.^2, dims=1)
	    # magnitudes squared of the tangential force
	    Ht_2 = sum(Ht.^2, dims=1)
	else # if w == Inf
		H0_x = [1.,0.,0.]
		psi_x = PotentialSimple_par(points, faces, normals, mu, H0_x)
	    Ht_x = HtField_par(points, faces, psi_x, normals)
	    Hn_x_norms = NormalFieldCurrent_par(points, faces, normals, Ht_x, mu, H0_x)
	    Hn_x = normals .* Hn_x_norms'

		H0_z = [0.,0.,1.]
		psi_z = PotentialSimple_par(points, faces, normals, mu, H0_z)
	    Ht_z = HtField_par(points, faces, psi_z, normals)
	    Hn_z_norms = NormalFieldCurrent_par(points, faces, normals, Ht_z, mu, H0_z)
	    Hn_z = normals .* Hn_z_norms'

	    # magnitudes squared of the normal force
	    Hn_2 = 0.5 * ( sum(Hn_x.^2, dims=1) + sum(Hn_z.^2, dims=1) )
	    # magnitudes squared of the tangential force
	    Ht_2 = 0.5 * ( sum(Ht_x.^2, dims=1) + sum(Ht_z.^2, dims=1) )
	end
	# f0 = mu*(mu-1)/8pi * Hn_2 + (mu-1)/8pi * Ht_2
	# f = f_0 + f_pert
	# efektīvi var teikt, ka:
	# Hn_2_pert = f_pert / ( mu*(mu-1)/8pi )
	if t == 0. #i == 1
		println("preturbing")
		alpha = 0.0 # perturbācijas izmērs
		k = 2 # cik staru

		phis = atan.(points[3,:],points[1,:])
		e_rhos = [[cos(phi),0.,sin(phi)] for phi in phis]
		fn_pert = [dot(e_rhos[k],normals[:,k]) for k in 1:size(points,2)]

		ampl = maximum(Hn_2)*mu*(mu-1)/8pi * alpha
		Hn_2_pert = -ampl * cos.(k*phis) .* fn_pert
		Hn_2 += Hn_2_pert'
	end


    @time velocities_phys = make_magvelocities_2_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
	@time velocities = make_Vvecs_conjgrad(normals,faces, points, velocities_phys, 1e-6, 500)

    #dt = 0.025*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
	if w != Inf
		dt = min(make_zinchencko_dt(points, connectivity, CDE, 7.4), (2*pi/w)*0.0505634375, 0.07)
	else
		dt = min(make_zinchencko_dt(points, connectivity, CDE, 7.4), 0.07)
	end
	# if i == 1 || (i - previous_i_when_split == 1 )
	# 	println("--- recalculating dt ---")
	# 	@time dt = make_meshdeg_dt(velocities, points, faces, connectivity; cutoff_crit = cutoff_crit - 0.01, degrade_after_steps = 5, eps = 10^-8)
	# 	if dt < 0
	# 		println("error: dt<0")
	# 		println("error: dt<0")
	# 		println("error: dt<0")
	# 		readline()
	# 	end
	# end
	println("---- dt = ", dt, " ----")
    t += dt
    points = points + velocities * dt

	#stdev = 0.001
	#r = randn(rng,Float64,size(points)) * stdev # add a gaussian noise with stdev
	#points += r

	if w != Inf
		H0 = [sin(w*t), 0., cos(w*t)]
	end

	# function vel_fun(points,t)
	# 	return make_magvelocities_2_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
	# end
	# function stab_fun(velocities_phys)
	# 	return make_Vvecs_conjgrad_par(normals,faces, points, velocities_phys, 1e-6, 500)
	# end

	#@time points, velocities_phys = rk2(points,t,dt,vel_fun, stab_fun; do_passive = true)


	# data = [points, faces]
	# println(points)
	# println(faces)
	# if i == 1
	# 	@save "hmm.jld2" data
	# end

    normals, CDE, AB = make_normals_parab(points, connectivity, normals; eps = 10^-8)

	# rescaling so that volume = V0
	rc = make_center_of_mass(points,faces, normals)
	V = make_volume(points,faces, normals)
	println("volume before rescale ", V)
	points = (points .- rc) * (V0/V)^(1/3) .+ rc
	println("rescaled to volume ", make_volume(points,faces, normals), "; start volume ", V0)
	println("center of mass ", rc)


    minN_triangles_to_split = 1

    marked_faces  = mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = cutoff_crit)
	println("Number of too large triangles: ",sum(marked_faces))
    if (sum(marked_faces) >= minN_triangles_to_split) && (i - previous_i_when_split >= 5)
        # split no more often than every 5 iterations to allow flipping to remedy the mesh
        previous_i_when_split = i

        # println("-----------------------------------")
        println("----------Adding mesh points-------")
        # println("    V-E+F = ", size(points,2)-size(edges,2)+size(faces,2))
        # println("    number of points: ", size(points,2))
        # println("    number of faces: ", size(faces,2))
        # println("    number of edges: ", size(edges,2))
        # println("-----------------------------------")

        global points_new, faces_new = add_points(points, faces,normals, edges, CDE; cutoff_crit = cutoff_crit)
        global edges_new = make_edges(faces_new)
        global connectivity_new = make_connectivity(edges_new)
		global marked_points = make_marked_points(points, faces, points_new, marked_faces)

        println("-----------------------------------")
        println("New V-E+F = ", size(points_new,2)-size(edges_new,2)+size(faces_new,2))
        println("New number of points: ", size(points_new,2))
        println("New number of faces: ", size(faces_new,2))
        println("New number of edges: ", size(edges_new,2))
        println("-----------------------------------")

		global do_active = true
		global oldE = 1000.
		print("locally relaxing added points")
		while do_active
		    points_new, maxgrad = relax_after_split_weights(marked_points, points_new, connectivity_new, points, normals, CDE,cutoff_crit; trianw = 5, trianpow = 2, edgepow = 2)

		    points_old = copy(points_new)

		    E = make_simple_E_weights(points_new,faces_new,points, normals, CDE,cutoff_crit; trianw = 5, trianpow = 2, edgepow = 2)
		    #println("Mesh energy ",E)
		    if E >= oldE
		        faces_new = faces_old
		        points_new = points_old
		        break
		    end

		    faces_old = copy(faces_new)
		    faces_new, connectivity_new, do_active = flip_edges(faces_new, connectivity_new, points_new)


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

    else # stabilize regularly if havent added new faces
        #H0 = [sin(w*t), 0., cos(w*t)]
        do_active = false

        faces, connectivity, do_active = flip_edges(faces, connectivity, points)

        if do_active
            if i - previous_i_when_flip > 5
                println("doing active because flipped")
                edges = make_edges(faces)
                points = active_stabilize(points, faces, CDE, connectivity, edges, normals)
            else
                println("flipped; not doing active")
            end
            previous_i_when_flip = i
        end

        if i % 100 == 0 && i > 2

            println("doing active every 100th time step")
            points = active_stabilize(points, faces, CDE, connectivity, edges, normals)

        end
    end

	a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])
    println(" --- c/a = $(c/a) , c/b = $(c/b)")

    if i % 1 == 0
		if w == Inf
			H0 = Hn = Ht = "not well defined"
		end
        data = [points, faces, t, velocities_phys, H0, Hn, Ht, normals, CDE, lambda, mu, w, Bm]
		#[points, faces, t, velocities_phys, H0, Bm]
        #println("Finished step $(last_step + i)")
        @save "$datadir/data$(lpad(i + last_step,5,"0")).jld2" data
    end
end # end simulation iterations
println("hooray, simulation finished :)")
