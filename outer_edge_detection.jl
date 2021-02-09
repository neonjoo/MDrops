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


dir = "star_5"
sourcedir = "/home/laigars/sim_data/$dir"

len = size(readdir(sourcedir),1) - 2

vs = []

dSs = []
Vs = zeros(Float64, len)

dira = readdir(sourcedir)[1500:1501]

len = size(dira, 1)
all_params = zeros(Float64, len, 10)
times = zeros(Float64, len)

global i, file
for (idx, file) in enumerate(dira)

    #dira = readdir(sourcedir)
    #last_file = readdir(sourcedir)[end-1200]
    println(idx, file)
    println("working on $file, i=$idx")
    @load "$sourcedir/$file" data
    points, faces, times[idx] = data[1], data[2], data[3]
    faces = Array{Int64,2}(faces)
    edges = make_edges(faces)
    connectivity = make_connectivity(edges)

    println(333)

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

        # if i == 464 || j == 464
        #     println()
        #     println("i=$i, j=$j, com=$common, prod=$(cross1[2]*cross2[2] > 0)")
        # end
        #
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
    #@save "outer2.jld2" outer_edges
    #%%
    #
    # function get_other_index(arr, ind)
    #     if arr[1] == ind
    #         return arr[2]
    #     else
    #         return arr[1]
    #     end
    # end
    # xs = zeros(size(uniques))
    # zs = zeros(size(uniques))
    #
    # 1151, 671
    # 671, 427
    # 427, 646   this ok
    # 646, 671
    # 671, 427

    #
    # pos = CartesianIndex(1,1)
    # for j in 1:size(uniques, 1)
    #     global pos
    #
    #     node = outer_edges[pos]
    #     other_node_pos = CartesianIndex(3 - pos[1], pos[2])
    #     other_node = outer_edges[other_node_pos]
    #
    #     println("$node, $other_node")#, other_node_pos)
    #     new_pos =  findall(x -> x == other_node, outer_edges)
    #     #println(new_pos)
    #     pos = get_other_index(new_pos, other_node_pos)
    #
    #     #new_edge = outer_edges[i]
    #     xs[j], zs[j] = points[1, node], points[3, node]
    #
    # end
    #%%
    ts = atan.(points[3, uniques], points[1,uniques])
    all = [ts points[1, uniques] points[3, uniques]]
    A = rand(1:100, 3, 4)   # a random matrix
    A[sortperm(A[:, 4]), :] # sorted by the 4th column

    alls = all[sortperm(all[:,1]), :]
    ts, xs, zs = alls[:, 1], alls[:, 2], alls[:, 3]
    rs = sqrt.(xs.^2 + zs.^2)

    ws = fft(rs)
    num_peaks = argmax(abs.(ws[2:20]))
    num_peaks = 6

    function polar(theta, params)
        R, a, w, phi = params
        return R .+ a .* cos.(w .* theta .+ phi)
    end

    function gaussian(theta, params)
        R, A, sigma, alpha = params[1:4]
        thetas = params[5:end]

        theta = repeat(theta, 1, size(thetas, 1))
        thetas = transpose(repeat(thetas, 1, size(theta, 1)))
        exp1 = thetas - theta
        exp2 = mod2pi.(thetas - theta .+ 2*pi)
        expon = abs.(cat(exp1, exp2, dims=3))
        expon = minimum(expon, dims=3)

        return vec(R .+ A .* sum(exp.(-1 ./sigma.^2 .* expon.^alpha), dims=2))
    end

    theta = range(0., stop=2*pi, length=200)
    #rr = polar(theta, 2, 1, 6, 0)

    #f(params) = sum((polar(ts, params) .- rs).^2)
    fg(params) = sum((gaussian(ts, params) .- rs).^2)

    #p0 = [1., 1., 6., 0.]
    phi = 0.
    peaks = range(0. + phi, 2*pi*(num_peaks-1)/num_peaks + phi, step=2*pi/num_peaks)
    p0g = vcat([minimum(rs), 0.5*(maximum(rs)-minimum(rs)), 0.3, 1.7], peaks)

    lower = vcat([0.8*minimum(rs), 0., 0.01, 1.], peaks .- 0.4)
    upper = vcat([maximum(rs), maximum(rs)-minimum(rs), 5., 4.], peaks .+ 0.4)

    #fits = optimize(f, p0)
    fitsg = optimize(fg, lower, upper, p0g)#, iterations=5000)
    params = fitsg.minimizer

    all_params[idx, :] = params

    scatter(ts, rs)
    plot!(ts, gaussian(ts, fitsg.minimizer), lw=2, show=true)
    #%%
    #scatter(rs .* cos.(ts), rs .* sin.(ts), lw=2)
    # plot!(gaussian(theta, p0g) .* cos.(theta), gaussian(theta, p0g) .* sin.(theta))
    #plot!(gaussian(theta, fitsg.minimizer) .* cos.(theta), gaussian(theta, fitsg.minimizer) .* sin.(theta), lw=3, show=true)

    #plot!(rr .* cos.(theta), rr .* sin.(theta))

    #%%
end

#%%
theta = range(0., stop=2*pi, length=200)
plot(1)
for p in all_params
    println(p)
    plot!(gaussian(theta, p) .* cos.(theta), gaussian(theta, p) .* sin.(theta), lw=2, show=true)
end
