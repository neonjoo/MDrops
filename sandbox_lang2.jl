using Plots
#using Makie
using JLD2
using FileIO

p = Plots

function eq431(e, mu)
    # 4pi ?
    N = 4pi*(1-e^2) / (2*e^3) * (log((1+e)/(1-e)) - 2e)
    temp1 = (3-2e^2)/e^2 - (3-4e^2)*asin(e)/(e^3 * sqrt(1-e^2))
    temp2 = (1-e^2)^(2/3) * ((3-e^2) * log((1+e)/(1-e))/e^5 - 6/e^4)

    # 4pi ?
    return (4pi/(mu-1) + N)^2 * 1/2pi * temp1/temp2
end

mu=3

dir = "elong_sphere_zinch4"
sourcedir = "/home/laigars/sim_data/$dir"
len = size(readdir(sourcedir),1) - 1

es = []
Bms = []
cs = []

for i in 5:10:50000
    # data = [points, faces, t, H0, Bm, v0max]

    @load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    #println("step $i")
    points = data[1]
    faces = data[2]
    Bm = data[end-1]

    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end

    #println(f([1.,1.]))
    x0 = [0.99, 1.01]
    res = Optim.optimize(f,x0)
    b = Optim.minimizer(res)[1]
    a = Optim.minimizer(res)[2]
    println("found $a and $b")
    
    c = (maximum(points[3,:]) - minimum(points[3,:]))/2
    a = (maximum(points[1,:]) + maximum(points[1,:]))/2

    e = sqrt(1-(a/c)^2)

    push!(es, e)
    push!(Bms, Bm)

    c = maximum(points[3,:])
    push!(cs, c)
end

p.scatter(es, 2pi*Bms)#, title="mu = $mu")
#title("mu = $mu")
p.xlabel!("e")
p.ylabel!("Bm")
p.plot!(es, eq431.(es, mu))
#p.plot!(es, eq431.(es, 30))
