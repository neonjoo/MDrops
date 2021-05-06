using JLD2
using FFTW
using Optim
#using LsqFit
using StatsBase
using LinearAlgebra

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")


dir = "perturbed_5_17_fastfield"

#sourcedir = "/home/laigars/sim_data/star_5"
sourcedir = "/mnt/hpc/sim_data/$dir"

len = size(readdir(sourcedir),1) - 2
#sourcedir = "/home/laigars/mdrops/meshes/perturbed/starfish_2"
dira = readdir(sourcedir)[3:2:end]

len = size(dira, 1)
println("size of $dira: $len")
all_params = zeros(Float64, len, 7)
times = zeros(Float64, len)

function gaussian(theta, params)
    #R, A, sigma, alpha = params[1:4]
    n = div(size(params, 1)-3, 2)
    R = params[1]
    A = params[2:n+1]
    sigma = params[n+2]
    alpha = params[n+3]
    thetas = params[n+4:end]

    #println(size(A,1))
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
wss = []
for (idx, file) in enumerate(dira)
    global wss
    #dira = readdir(sourcedir)
    #last_file = readdir(sourcedir)[end-1200]
    println(idx, file)
    println("working on $file, i=$idx")
    @load "$sourcedir/$file" data
    #@load "$file" data
    points, faces = data[1], data[2]
    #times[idx] = data[3]
    faces = Array{Int64,2}(faces)
    edges = make_edges(faces)
    connectivity = make_connectivity(edges)

    println(333)
    global pars
    #%%
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

    global ws = fft(rs)
    push!(wss, abs.(ws))
    num_peaks = argmax(abs.(ws[2:20]))
    #num_peaks = 6

    function polar(theta, params)
        R, a, w, phi = params
        return R .+ a .* cos.(w .* theta .+ phi)
    end


    theta = range(0., stop=2*pi, length=200)
    #rr = polar(theta, 2, 1, 6, 0)

    #f(params) = sum((polar(ts, params) .- rs).^2)
    fg(params) = sum((gaussian(ts, params) .- rs).^2)


    #p0 = [1., 1., 6., 0.]
    phi = 0.
    peaks = range(0. + phi, 2*pi*(num_peaks-1)/num_peaks + phi, step=2*pi/num_peaks)


    global p0g = vcat([0.8*minimum(rs)], repeat([0.9*(maximum(rs)-minimum(rs))], num_peaks), [1., 1.9], peaks)


    if idx > 1
        p0g = params
    end

    lower = vcat([0.5*minimum(rs)], repeat([0.], num_peaks), [0.01, 1.], peaks .- 0.4)
    upper = vcat([1.03*minimum(rs)], repeat([2*(maximum(rs)-minimum(rs))], num_peaks), [5., 4.], peaks .+ 0.4)

    #fits = optimize(f, p0)
    println("optimizing")
    #@time fitsg = optimize(fg, lower, upper, p0g)#, iterations=5000)
    #global params = fitsg.minimizer

    #all_params[idx, :] = params

    Plots.scatter(ts, rs)
    #plot!(ts, gaussian(ts, fitsg.minimizer), lw=2, show=true)
    #%%
    #scatter(rs .* cos.(ts), rs .* sin.(ts), lw=2)
    # plot!(gaussian(theta, p0g) .* cos.(theta), gaussian(theta, p0g) .* sin.(theta))
    #plot!(gaussian(theta, fitsg.minimizer) .* cos.(theta), gaussian(theta, fitsg.minimizer) .* sin.(theta), lw=3, show=true)

    #plot!(rr .* cos.(theta), rr .* sin.(theta))

    #%%
end

#%%


#theta = range(0., stop=2*pi, length=200)

#for p in all_params
# #    println(p
# scatter(xs, zs)
# plot!(gaussian(theta, all_params[end, :]) .* cos.(theta), gaussian(theta, all_params[end, :]) .* sin.(theta), lw=2, show=true)
# plot(ts, rs)
# plot!(ts, gaussian(ts, all_params[end, :]), lw=2, show=true)
#
#
# all_times = [times17, times19, times21, times23, times30]
# all_p = [p17, p19, p21, p23, p30]
#
# data = [all_times, all_p]
# @save "./perturbation_times.jld2" all_times
# @save "./perturbation_params.jld2" all_p

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
