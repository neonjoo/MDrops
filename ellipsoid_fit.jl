using Makie
using JLD2
using FileIO
using Optim

dir = "2019-03-26/3"
#sourcedir = "/home/andris/sim_data/$dir"
sourcedir = "/home/joo/sim_data/collapse_ellipse_zinch8/"

len = size(readdir(sourcedir),1) - 1


ts = zeros(Float64,len)
a_over_b = zeros(Float64,len)
D = zeros(Float64,len)
for i in 1:len
    global ts, a_over_b
    @load "$sourcedir/data$(lpad(i,5,"0")).jld2" data

    points = data[1]
    faces = data[2]
    t = data[3]

    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end

    #println(f([1.,1.]))
    x0 = [0.99, 1.01]
    res = Optim.optimize(f,x0)
    b = Optim.minimizer(res)[1]
    a = Optim.minimizer(res)[2]

    ts[i] = t
    a_over_b[i] = a/b
    D[i] = (a-b)/(a+b)

end

using PyPlot
PyPlot.plot(ts[1:end],D[1:end],".")

dir = "2019-03-26/2"
#sourcedir = "/home/andris/sim_data/$dir"
len2 = size(readdir(sourcedir),1) - 1
#len=0

#ts = append!(ts,zeros(Float64,len2))
#a_over_b = append!(a_over_b,zeros(Float64,len2))
#D = append!(D,zeros(Float64,len2))
# ts = zeros(Float64,len2)
# a_over_b = zeros(Float64,len2)
for i in 1:len2

    global ts, a_over_b
    @load "$sourcedir/data$(lpad(i,5,"0")).jld2" data

    points = data[1]
    faces = data[2]
    t = data[3]

    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end

    println(f([1.,1.]))
    x0 = [0.99, 1.01]
    res = Optim.optimize(f,x0)
    b = Optim.minimizer(res)[1]
    a = Optim.minimizer(res)[2]

    ts[i] = t
    a_over_b[i] = a/b
    D[i] = (a-b)/(a+b)

end

using PyPlot

pygui()
figure()

p.plot(ts[1:end],D[1:end])#,".")
lambda = 1.
τ = (16+19*lambda)*(3+2*lambda)/(40+40*lambda)
logD0 = log(D[end]) + ts[end]/τ
teorlogD = logD0 .- ts/τ

p.plot(ts[1:end],teorlogD[1:end], label="teor")#,"-")
p.title!("$sourcedir")
p.plot!(ts[1:end],log.(D[1:end]), label="num")
#p.plot(vec(1:length(sqrt.(sum(velocities.^2,dims=1)))),sort(vec(sqrt.(sum(velocities.^2,dims=1)))))
