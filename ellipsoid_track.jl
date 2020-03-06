using Makie
using JLD2
using FileIO
using Optim
using LinearAlgebra
using StatsBase

#sourcedir = "/home/andris/sim_data/$dir"
sourcedir = "/home/laigars/sim_data/rotating_2.3"

len = size(readdir(sourcedir),1) - 2

println("Number of files in $sourcedir: $len")
println()

ts = zeros(Float64, len)
abcs = zeros(3, len)
eigvecs_z = zeros(3, len)
angles_z = zeros(1, len)

for i in 1:50:len
    global ts, points, faces
    @load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    #@load "./meshes/faces_critical_hyst_2_21.jld2" faces
    #@load "./meshes/points_critical_hyst_2_21.jld2" points
    #points = Array{Float64}(transpose(points))
    points = data[1]
    faces = data[2]
    t = data[3]
    print("$i ")
    #mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
    #                StatsBase.mean(points[2,:]),
    #                StatsBase.mean(points[3,:]))

    #points = points .- [mean_x, mean_y, mean_z]

    function f(coefs::Array{Float64,1})
        # coefs - alphas and betas for athe ellipsoid equation
        a1, a2, a3, b1, b2, b3 = coefs
        x,  y, z = points[1,:], points[2,:], points[3,:]
        eq = a1*x.^2 + a2*y.^2 + a3*z.^2 + 2*b1*y.*z + 2*b2*x.*z + 2*b3*y.*z .- 1

        return sum(eq.^2)
    end

    x0 = [0.99, 1.01, 1., 0., 0., 0.]
    res = Optim.optimize(f,x0)
    a1, a2, a3, b1, b2, b3 = Optim.minimizer(res)
    global M = [a1 b3 b2; b3 a2 b1; b2 b1 a3] # of the X^T * M * X = 1 equation for a transformed ellipsoid
    #println("koefs: ", Optim.minimizer(res))
    #println()
    #println("vecs: ", transpose(eigvecs(M)))
    if i == 87
        println()
        scene = Makie.mesh(points', faces',color = :gray, shading = false, visible = true)
        Makie.wireframe!(scene[end][1], color = :black, linewidth = 0.7,visible = true)
    end

    eigs = eigvals(M)
    vecs = eigvecs(M)
    eigvecs_z[:, i] = vecs[:, 1]
    angles_z[i] = acos(dot(vecs[:,1], [0., 0., 1.])) / pi * 180
    # according to
    # https://math.stackexchange.com/questions/944336/equation-of-major-axis-of-an-ellipsoid
    abcs[:, i] = sqrt.(1 ./ abs.(eigs))

end

#
i = 301
@load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
points = data[1]
faces = data[2]

scene = Makie.mesh(points', faces',color = :gray, shading = false, visible = true)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 0.7,visible = true)#, limits=FRect3D((-5,-5,-5),(10,10,10)))
println(abcs[:, i])

p.plot(angles_z')
