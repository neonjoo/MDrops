using Pkg
pkg"activate ."
pkg"resolve"

using Plots
#using Makie
using JLD2
using FileIO
using Optim
using CSV

p = Plots
p.pyplot()

include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry

function eq431(e, mu)
    # 4pi ?
    e = sqrt(1-e^-2)
    N = 4pi*(1-e^2) / (2*e^3) * (log((1+e)/(1-e)) - 2e)
    temp1 = (3-2e^2)/e^2 - (3-4e^2)*asin(e)/(e^3 * sqrt(1-e^2))
    temp2 = (1-e^2)^(2/3) * ((3-e^2) * log((1+e)/(1-e))/e^5 - 6/e^4)

    # 4pi ?
    return (4pi/(mu-1) + N)^2 * 1/2pi * temp1/temp2
end

function shank(seq)
    return seq[3] - (seq[3] - seq[2])^2 / ((seq[3] - seq[2]) - (seq[2] - seq[1]))
end

mu=35.5

dir = "elong_sphere_12"
sourcedir = "/home/laigars/sim_data/$dir"
len = size(readdir(sourcedir),1) - 1

ba_exp = CSV.read("/home/laigars/sim_data/bas.csv", header=0)[1]
t_exp = CSV.read("/home/laigars/sim_data/time.csv", header=0)[1]

es = []
num_es = []
Bms = []
num_Bms = []
cs = []
all_Bms = []
all_es = []
last_Bm = 0
last_t = 0

volumes = []
ts = []
all_ts = []
for f in readdir(sourcedir)[2:end]
    # data = [points, faces, t, H0, Bm, v0max]
    global last_Bm, es
    #@load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    @load "$sourcedir/$f" data
    #println("step $i")
    points = data[1]
    faces = data[2]
    Bm = data[end-1]
    t = data[3]
    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end

    #println(f([1.,1.]))
    x0 = [0.99, 1.01]
    res = Optim.optimize(f,x0)
    a = Optim.minimizer(res)[1]
    b = Optim.minimizer(res)[2]

    e = sqrt(1-(a/b)^2)
    e = b/a

    push!(num_es, e)
    #push!(num_Bms, Bm)
    push!(volumes, -SG.volume(points, faces))
    push!(ts, t)
    #
    # if Bm > last_Bm
    #     es = []
    #     push!(all_Bms, Bm)
    #     last_Bm = Bm
    #     push!(all_es, es)
    #     push!(all_ts, t)
    # end
    push!(es, e)

end
push!(all_ts, ts[end])


#p.scatter(ts .+ t_exp[34], num_es, label="parameters from fit 1", color="red",markersize=1,markerstrokecolor="red", show=true)


dir = "elong_sphere_14"
sourcedir = "/home/laigars/sim_data/$dir"
len = size(readdir(sourcedir),1) - 1

ba_exp = CSV.read("/home/laigars/sim_data/bas.csv", header=0)[1]
t_exp = CSV.read("/home/laigars/sim_data/time.csv", header=0)[1]

es = []
num_es = []
Bms = []
num_Bms = []
cs = []
all_Bms = []
all_es = []
last_Bm = 0
last_t = 0

volumes = []
ts = []
all_ts = []
for f in readdir(sourcedir)[2:end]
    # data = [points, faces, t, H0, Bm, v0max]
    global last_Bm, es
    #@load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    @load "$sourcedir/$f" data
    #println("step $i")
    points = data[1]
    faces = data[2]
    Bm = data[end-1]
    t = data[3]
    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end

    #println(f([1.,1.]))
    x0 = [0.99, 1.01]
    res = Optim.optimize(f,x0)
    a = Optim.minimizer(res)[1]
    b = Optim.minimizer(res)[2]

    e = sqrt(1-(a/b)^2)
    e = b/a

    push!(num_es, e)
    #push!(num_Bms, Bm)
    push!(volumes, -SG.volume(points, faces))
    push!(ts, t)
    #
    # if Bm > last_Bm
    #     es = []
    #     push!(all_Bms, Bm)
    #     last_Bm = Bm
    #     push!(all_es, es)
    #     push!(all_ts, t)
    # end
    push!(es, e)

end


p.scatter(ts .+ t_exp[34], num_es, label="lower bound", color="green",markersize=1,markerstrokecolor="green", show=true)


dir = "elong_sphere_13"
sourcedir = "/home/laigars/sim_data/$dir"
len = size(readdir(sourcedir),1) - 1

ba_exp = CSV.read("/home/laigars/sim_data/bas.csv", header=0)[1]
t_exp = CSV.read("/home/laigars/sim_data/time.csv", header=0)[1]

es = []
num_es = []
Bms = []
num_Bms = []
cs = []
all_Bms = []
all_es = []
last_Bm = 0
last_t = 0

volumes = []
ts = []
all_ts = []
for f in readdir(sourcedir)[2:end]
    # data = [points, faces, t, H0, Bm, v0max]
    global last_Bm, es
    #@load "$sourcedir/data$(lpad(i,5,"0")).jld2" data
    @load "$sourcedir/$f" data
    #println("step $i")
    points = data[1]
    faces = data[2]
    Bm = data[end-1]
    t = data[3]
    function f(ab::Array{Float64,1})
        return sum((points[1,:].^2/ab[1]^2 .+ points[2,:].^2/ab[1]^2 .+
                points[3,:].^2/ab[2]^2 .- 1).^2)
    end

    #println(f([1.,1.]))
    x0 = [0.99, 1.01]
    res = Optim.optimize(f,x0)
    a = Optim.minimizer(res)[1]
    b = Optim.minimizer(res)[2]

    e = sqrt(1-(a/b)^2)
    e = b/a

    push!(num_es, e)
    #push!(num_Bms, Bm)
    push!(volumes, -SG.volume(points, faces))
    push!(ts, t)
    #
    # if Bm > last_Bm
    #     es = []
    #     push!(all_Bms, Bm)
    #     last_Bm = Bm
    #     push!(all_es, es)
    #     push!(all_ts, t)
    # end
    push!(es, e)

end





push!(all_ts, ts[end])
# finals = []
# for arr in num_es
#     final = shank(arr[end-2:end])
#     push!(finals, final)
# end
finals = []

for i in 1:37
    final = shank(num_es[end-(1+i):end-(i-1)])
    push!(finals, final)
end

p.scatter!(ts .+ t_exp[34], num_es, label="upper bound", color="blue",markersize=1,markerstrokecolor="blue", show=true)
p.plot!(t_exp, ba_exp, label="experimental", color="black", linewidth=2, framestyle=:box)
#p.plot!([0,30], [sb.mean(finals), sb.mean(finals)], label="avg shanked")

#p.plot!([0,30], [finals[1], finals[1]], label="last shanked")


# p.scatter(num_es, num_Bms, label="all_many", color="lightblue")
# p.scatter!(finals, all_Bms, label="shanked_many", color="blue")#, title="mu = $mu")

#p.scatter!(num_es_few, num_Bms_few, label="all_few", color="orange")
#p.scatter!(finals_few, all_Bms_few, label="shanked_few", color="red")#, title="mu = $mu")
#title("mu = $mu")
p.xlabel!("t")
p.ylabel!("b/a")
p.xlims!(0,140)
p.ylims!(1,3)
p.plot!([25,25], [1,2.8],color="black", linestyle=:dash, label="", legend=:topleft)
p.annotate!(14, 1.17, text("1.178 Oe", 8))

p.plot!([55,55], [1,2.8],color="black", linestyle=:dash, label="")
p.annotate!(45, 1.87, text("1.872 Oe", 8))

p.plot!([84, 84], [1,2.8],color="black", linestyle=:dash, label="")
p.annotate!(74, 2.3, text("2.344 Oe", 8))

#e = 0:0.002:finals[end]+0.01p
p.savefig("num_vs_exp.pdf")

e=0.01:0.001:0.90

ba = 1:0.01:30
#p.plot(ba, eq431.(ba, mu), legend=:right, label="teor", color="green")

#p.plot(e, (3/4/pi*4.141)^(1/3) * eq431.(e, 3), legend=:right, title=:"mu=3", label="teor_scaled")
#p.plot!(es, eq431.(es, 30))
