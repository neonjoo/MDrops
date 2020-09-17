using Pkg
pkg"activate ."
pkg"resolve"

#using Makie
using FileIO
using JLD2
using StatsBase
using Optim
#
# datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30/"
#
# files = readdir(datadir)
#
#
#
# Ndata = size(files,1)-3
# as = zeros(Ndata,1)
# bs = zeros(Ndata,1)
# ts = zeros(Ndata,1)
# for i = 1:Ndata
#     file = files[2+i]
#     #println(file)
#     @load "$datadir/$file" data
#
#     points, faces, t = data[1], data[2], data[3]
#     faces = Array{Int64,2}(faces)
#     mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
#                     StatsBase.mean(points[2,:]),
#                     StatsBase.mean(points[3,:]))
#     points = points .- [mean_x, mean_y, mean_z]
#
#     function f(ab::Array{Float64,1})
#         return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
#                 points[3,:].^2/ab[2]^2 .- 1).^2)
#     end
#     x0 = [0.99, 1.01]
#     res = Optim.optimize(f,x0)
#     b = Optim.minimizer(res)[1]
#     a = Optim.minimizer(res)[2]
#
#     as[i] = a
#     bs[i] = b
#     ts[i] = t
# end

using Plots
#plot(ts,as./bs,legend=:bottom)

datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_adaptive_dt/"

files = readdir(datadir)

Ndata = size(files,1)-3
#Ndata = 20
as = zeros(Ndata,1)
bs = zeros(Ndata,1)
ts = zeros(Ndata,1)
for i = 1:Ndata
    file = files[2+i]
    #println(file)
    @load "$datadir/$file" data

    points, faces, t = data[1], data[2], data[3]
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

    as[i] = a
    bs[i] = b
    ts[i] = t
end

plot(ts,as./bs)

datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_manyN_adaptive_dt/"

files = readdir(datadir)

Ndata = size(files,1)-3
#Ndata = 50
as = zeros(Ndata,1)
bs = zeros(Ndata,1)
ts = zeros(Ndata,1)
for i = 1:Ndata
    file = files[2+i]
    #println(file)
    @load "$datadir/$file" data

    points, faces, t = data[1], data[2], data[3]
    faces = Array{Int64,2}(faces)
    mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                    StatsBase.mean(points[2,:]),
                    StatsBase.mean(points[3,:]))
    points = points .- [mean_x, mean_y, mean_z]

    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end
    x0 = [0.5, 1.5]
    res = Optim.optimize(f,x0)
    b = Optim.minimizer(res)[1]
    a = Optim.minimizer(res)[2]

    as[i] = a
    bs[i] = b
    ts[i] = t
end

plot!(ts,as./bs)


datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_manymoreN_adaptive_dt/"

files = readdir(datadir)

Ndata = size(files,1)-3
#Ndata = 50
as = zeros(Ndata,1)
bs = zeros(Ndata,1)
ts = zeros(Ndata,1)
for i = 1:Ndata
    file = files[2+i]
    #println(file)
    @load "$datadir/$file" data

    points, faces, t = data[1], data[2], data[3]
    faces = Array{Int64,2}(faces)
    mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                    StatsBase.mean(points[2,:]),
                    StatsBase.mean(points[3,:]))
    points = points .- [mean_x, mean_y, mean_z]

    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end
    x0 = [0.5, 1.5]
    res = Optim.optimize(f,x0)
    b = Optim.minimizer(res)[1]
    a = Optim.minimizer(res)[2]

    as[i] = a
    bs[i] = b
    ts[i] = t
end


plot!(ts,as./bs)



datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_adaptiveN_adaptive_dt/"

files = readdir(datadir)

Ndata = size(files,1)-3
#Ndata = 50
as = zeros(Ndata,1)
bs = zeros(Ndata,1)
ts = zeros(Ndata,1)
for i = 1:Ndata
    file = files[2+i]
    #println(file)
    @load "$datadir/$file" data

    points, faces, t = data[1], data[2], data[3]
    faces = Array{Int64,2}(faces)
    mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                    StatsBase.mean(points[2,:]),
                    StatsBase.mean(points[3,:]))
    points = points .- [mean_x, mean_y, mean_z]

    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end
    x0 = [0.5, 1.5]
    res = Optim.optimize(f,x0)
    b = Optim.minimizer(res)[1]
    a = Optim.minimizer(res)[2]

    as[i] = a
    bs[i] = b
    ts[i] = t
end


plot!(ts,as./bs)
