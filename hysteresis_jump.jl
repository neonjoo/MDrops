using Pkg
pkg"activate ."
pkg"resolve"

using Plots
#using Makie
using JLD2
using FileIO
using Optim
using StatsBase
p = Plots

include("./SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry

function eq431_ab(e, mu)
    # 4pi ?
    e = sqrt(1-e^-2)
    N = 4pi*(1-e^2) / (2*e^3) * (log((1+e)/(1-e)) - 2e)
    temp1 = (3-2e^2)/e^2 - (3-4e^2)*asin(e)/(e^3 * sqrt(1-e^2))
    temp2 = (1-e^2)^(2/3) * ((3-e^2) * log((1+e)/(1-e))/e^5 - 6/e^4)

    # 4pi ?
    return (4pi/(mu-1) + N)^2 * 1/2pi * temp1/temp2
end

function eq431_e(e, mu)
    # 4pi ?
    #e = sqrt(1-e^-2)
    N = 4pi*(1-e^2) / (2*e^3) * (log((1+e)/(1-e)) - 2e)
    temp1 = (3-2e^2)/e^2 - (3-4e^2)*asin(e)/(e^3 * sqrt(1-e^2))
    temp2 = (1-e^2)^(2/3) * ((3-e^2) * log((1+e)/(1-e))/e^5 - 6/e^4)

    # 4pi ?
    return (4pi/(mu-1) + N)^2 * 1/2pi * temp1/temp2
end


function shank(seq)
    return seq[3] - (seq[3] - seq[2])^2 / ((seq[3] - seq[2]) - (seq[2] - seq[1]))
end

mu=30

dir = "hysteresis_1_return3"
sourcedir = "/home/laigars/sim_data/$dir"
len = size(readdir(sourcedir),1) - 1


es = []
num_es = []
num_bas = []
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

mean_ps = zeros(3, len)
i=1
for f in readdir(sourcedir)[2:20:end-1]
    global last_Bm, es, i
    @load "$sourcedir/$f" data
    points = data[1]
    faces = data[2]
    Bm = data[end-1]
    t = data[3]

    mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                    StatsBase.mean(points[2,:]),
                    StatsBase.mean(points[3,:]))

    mean_ps[:,i] = [mean_x, mean_y, mean_z]
    points = points .- mean_ps[:,i]

    i += 1
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
    ba = b/a

    push!(num_es, e)
    push!(num_bas, ba)
    push!(num_Bms, Bm)
    push!(volumes, -SG.volume(points, faces))
    push!(ts, t)
    if Bm > last_Bm
        es = []
        push!(all_Bms, Bm)
        last_Bm = Bm
        push!(all_es, es)
        push!(all_ts, t)
    end
    push!(es, e)

end

finals = []
for arr in all_es
    final = shank(arr[end-2:end])
    push!(finals, final)
end
finals = 1 ./ (1 .- sqrt.(finals))

#p.plot!([0,30], [sb.mean(finals), sb.mean(finals)], label="avg shanked")

num_es = Array{Float64,1}(num_es)
num_bas = Array{Float64,1}(num_bas)
num_Bms = Array{Float64,1}(num_Bms)
finals = Array{Float64,1}(finals)
all_Bms = Array{Float64,1}(all_Bms)
ts = Array{Float64,1}(ts)
es = Array{Float64,1}(es)

p.scatter(num_bas, num_Bms, label="numerical", color="blue", markerstrokecolor=:lightblue)
#p.scatter!(finals, all_Bms, label="shanked_many", color="red")#, title="mu = $mu")
#p.plot(ts, num_es, label="")
#p.plot!(ts, num_Bms./5, label="")
#p.scatter!(num_es_few, num_Bms_few, label="all_few", color="orange")
#p.scatter!(finals_few, all_Bms_few, label="shanked_few", color="red")#, title="mu = $mu")
#title("mu = $mu")
p.plot(num_bas)
p.xlabel!("b/a")
p.ylabel!("Bm")
#p.xlims!(1,12)
p.ylims!(4,7)
p.title!("$sourcedir")
#e = 0:0.002:finals[end]+0.01
#p.xticks!(1:0.5:10)
p.yticks!(1:0.1:9)
bas = [7.6, 9.8, 1.85, 1.36, 1.17, 1.06]
bas2 = [8.01, 7.03, 5.69, 5.175]

bmss = [4, 5, 3.8, 2.8, 1.8, 0.8]
bmss2 = [3.8, 3.4, 3.2, 3.0]

p.scatter(bas, bmss, color=:blue, markersize=8, label="increasing Bm")
p.scatter!(bas2, bmss2, color=:red, markersize=8, label="decreasing Bm")
e=0.01:0.001:0.999
ba = 1:0.01:10
#p.plot(ba, eq431.(ba, mu), legend=:top, label="teor", color="green")
#p.plot(e, (3/4/pi*4.141)^(1/3) * eq431.(e, 30), legend=:right, title=:"mu=3", label="teor_scaled")
p.plot!(ba, eq431_ab.(ba, 30), color=:green, label="theoretical")







e_c = 0.8926
ba_c = 2.2189
Bm_crit = 3.68423
es_crit = num_es[num_es .> e_c]
ts_crit = ts[end-size(es_crit,1)+1:end]

ind0 = findfirst(x -> x > 1100, ts)
ind1 = findfirst(x -> x > 8, ts)
ba0 = findfirst(x -> x > ba_c, num_bas)

ba_fit = num_bas[ba0:ind1]
ts_fit = ts[ba0:ind1]
function tange(pp::Array{Float64,1})
    return sum((ba_fit .- ba_c .- pp[1]*pp[2] * tan.((ts_fit .- ts_fit[1])./pp[2])).^2)
end


x0 = [0.0001, 100]
res = Optim.optimize(tange,x0)
S = Optim.minimizer(res)[1]
tau = Optim.minimizer(res)[2]

#p.plot(ts_crit, es_crit, label="num > e_c (= 0.8923)", legend=:left, color=:black)
p.plot!(ts, num_bas, label="num_all", legend=:right)
p.title!("mu=30")
#p.plot!([1, 30], [ba_c, ba_c], label="ab_critical (2.219)")
p.xlabel!("time, s")
p.ylabel!("b/a")
#p.plot!(ts_s, es_s, label="exp_fit")
#p.ylims!(2, 4)
#p.xlims!(0, 30)
p.plot!(ts_fit, ba_c .+ S * tau * tan.((ts_fit .- ts_fit[1])./tau), label="fit")

#p.plot!(ts_fit, ba_c .+ S * 10*tau * tan.((ts_fit .- ts_fit[1])./tau/10), label="fit2")
p.title!("S = $S, tau = $tau")

println(tan((ts_fit[end] - ts_fit[1])/tau))
