using Makie
using JLD2
using FileIO
using Optim

dir = "pushing_to_limit_langfix_extended"
sourcedir = "/home/andris/sim_data/$dir"
len = size(readdir(sourcedir),1) - 1


ts = zeros(Float64,len)
a_over_b = zeros(Float64,len)
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

end

using PyPlot
PyPlot.plot(ts,a_over_b,".")
