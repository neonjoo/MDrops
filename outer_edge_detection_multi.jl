using JLD2
using FFTW
using Optim
#using LsqFit
using StatsBase
using LinearAlgebra
using Distributed
using Plots

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")

bms = [17, 18, 19, 20, 21, 22, 23, 24, 25]
#bms = [21, 22]
#bms = [21]
ns_pars = []

for n in 6:1:6

    global all_pars = []
    global all_ts = []
    global all_ws = []
    for bm in bms

        outer_dir = "/mnt/hpc/sim_data"#perturbed_$(n)"
        dir = "perturbed_$(n)_$(bm)_fastfield"
        sourcedir = "$outer_dir/$dir"
        println("-------------------------------------------------------------------- $dir")

        global wss = []

        len = size(readdir(sourcedir),1) - 3
        start = 3
        step = 1
        dira = readdir(sourcedir)[start:step:end]
        println(dira)

        all_params = zeros(Float64, size(dira,1), 3 + 2*n)
        times = zeros(Float64, size(dira,1))

        function gaussian(theta, params)
            #R, A, sigma, alpha = params[1:4]
            n = div(size(params, 1)-3, 2)
            R = params[1]
            A = params[2:n+1]
            sigma = params[n+2]
            alpha = params[n+3]
            thetas = params[n+4:end]

            #println("size of A = $(size(A,1))")
            #println(size(theta,1))
            As = transpose(repeat(A, 1, size(theta, 1)))
            theta = repeat(theta, 1, size(thetas, 1))
            thetas = transpose(repeat(thetas, 1, size(theta, 1)))
            exp1 = thetas - theta
            exp2 = mod2pi.(thetas - theta .+ 2*pi)
            expon = abs.(cat(exp1, exp2, dims=3))
            expon = minimum(expon, dims=3)
            gaus = exp.(-1 ./sigma.^2 .* expon.^alpha)
            #display(gaus)

            #display(As)

            return vec(R .+ sum(As .* gaus, dims=2))
            #return vec(R .+ sum(A .* exp.(-1 ./sigma.^2 .* expon.^alpha), dims=2))
        end

        global i, file
        for (idx, file) in enumerate(dira)

            println("working on $file, Bm=$bm, i=$idx")
            @load "$sourcedir/$file" data
            #@load "$file" data
            points, faces = data[1], data[2]

            times[idx] = data[3]
            faces = Array{Int64,2}(faces)

            if idx == 1
                print("making first edges.. -----------")
                @time global edges = make_edges(faces)
                global Npoints = size(points, 2)
            else
                if size(points,2) != Npoints
                    print("remaking edges.. ################")
                    @time global edges = make_edges(faces)
                    Npoints = size(points, 2)
                end # end if new_points = old_points
                #println("same points, moving on at step $idx !")
            end # end if idx=1

            print("making connectivity.. ")
            @time global connectivity = make_connectivity(edges)

            function is_outer(i, j, points, connectivity)
                # Return if edge between nodes i--j is an outer edge

                # project points into XZ plane
                points[2,:] .= 0

                # find the 2 common neighbors n1 and n2
                common = intersect(connectivity[:,i], connectivity[:,j])
                common = filter(x -> x>0, common)

                r_ij = points[:,j] - points[:,i]
                r_in1 = points[:,common[1]] - points[:,i]
                r_in2 = points[:, common[2]] - points[:,i]

                cross1 = cross(r_ij, r_in1)
                cross2 = cross(r_ij, r_in2)

                return cross1[2]*cross2[2] > 0
            end

            global outer_edges = [0;0]

            for k in 1:size(edges, 2)
                global outer_edges
                i, j = edges[1,k], edges[2,k]
                if is_outer(i, j, points, connectivity)
                    outer_edges = hcat(outer_edges, [i; j])
                end
            end

            outer_edges = outer_edges[:,2:end]
            uniques = unique(outer_edges) # the unique nodes

            global ts = atan.(points[3, uniques], points[1,uniques]) .+ pi
            all = [ts points[1, uniques] points[3, uniques]]

            global alls = all[sortperm(all[:,1]), :]
            global ts, xs, zs, ws
            ts, xs, zs = alls[:, 1], alls[:, 2], alls[:, 3]
            global rs = sqrt.(xs.^2 + zs.^2)

            ws = fft(rs)
            num_peaks = argmax(abs.(ws[2:10]))
            push!(wss, abs.(ws))
            println(idx, wss[end][3])
            theta = range(0., stop=2*pi, length=200)
            fg(params) = sum((gaussian(ts, params) .- rs).^2)

            phi = 0.
            peaks = range(0. + phi, 2*pi*(num_peaks-1)/num_peaks + phi, step=2*pi/num_peaks)


            global p0g = vcat([0.8*minimum(rs)], repeat([0.9*(maximum(rs)-minimum(rs))], num_peaks), [1., 1.9], peaks)

            if idx > 1
                p0g = params .* 0.85
            end

            global lower = vcat([0.5*minimum(rs)], repeat([0.], num_peaks), [0.01, 1.], peaks .- 0.4) .* 0.8
            global upper = vcat([1.03*minimum(rs)], repeat([2*(maximum(rs)-minimum(rs))], num_peaks), [5., 4.], peaks .+ 0.4) .* 1.5
            print("optimizing: ")
            #@time fitsg = optimize(fg, lower, upper, p0g)#, iterations=5000)
            #display(params)
            #println()
            #global params = fitsg.minimizer

            #all_params[idx, :] = params

            # Plots.scatter(ts, rs, title="$file, Bm=$bm, i=$idx")
            # plot!(ts, gaussian(ts, params), lw=2)#, show=true)
            # savefig("~/sim_data/perturbations/$(n)_$(bm)_$file.png")

        end # end folder loop

        push!(all_ts, times)
        push!(all_ws, wss)
        push!(all_pars, all_params)

        @save "./perturbations/perturbation_fourier_n_$(n)_Bm_$(bm).jld2" wss
    end # end bm loop

    #@save "./perturbation_times_fastfield_n_$n.jld2" all_ts
    #@save "./perturbation_params_n_$n.jld2" all_pars
    #@save "./perturbation_fourier_fastfield_n_$n.jld2" all_ws
end  # end n loop
#%%


#
# data = [all_times, all_p]


#%%
# fs = 12
#
# plot(abs.(ws[2:10]), label="", lw=2, color=1)
# Plots.scatter!(abs.(ws[2:10]), label="", color=1, markersize=5)
# Plots.ylabel!("A", labelfontsize=fs)
# Plots.xlabel!("n", labelfontsize=fs)
#
# #%%
#
# params = [2.4293003, 0.4943, 0.86815, 0.9672, 0.9611, 0.9821, 1.185, 0.2435, 1.7895, 0.2159, 1.203, 2.1892, 3.3028, 4.2414, 5.408]
#
# params = all_params[end,:]
#
# #all_params = params
# scatter(xs, zs, label="")
#
# #theta = theta .- pi
# plot!(gaussian(theta, params) .* cos.(theta), gaussian(theta, params) .* sin.(theta), lw=2, show=true, color=:black, label="")
# #
# scatter(ts, rs)
# plot!(ts, gaussian(ts, all_params[end, :]), lw=2, show=true)
# #%%
#
# pars = [1.221, 0.039, 0.054, 0.038, 0.051, 0.051, 0.260, 2.214, 0.016, 1.263, 2.489, 3.848, 5.059]
# sh = 0

# plot(ts .- sh, gaussian(ts .- sh, pars), lw=2, show=true)

#%%

# plot()
# for (i, ww) in enumerate(all_ws)
#     global qq = []
#     for w in ww[2:end]
#         global qq
#         push!(qq, w[3])
#
#     end
#     qq = Array{Float64}(qq)
#     println(qq)
#
#     p = plot!(all_ts[i][2:end-2], qq, label=bms[i])
#     display(p)
# end
