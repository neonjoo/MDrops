using Pkg
pkg"activate ."
pkg"resolve"

using Plots
#using Makie
using JLD2
using FileIO
using Optim
using CSV

#p = Plots
#p.pyplot()
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")

points, faces = expand_icosamesh(R=1, depth=2)
points = Array{Float64}(points)
faces = Array{Int64}(faces)
normals = Normals(points, faces)

function make_deltaH_normal(points, faces, alpha, H0)

    normals = Normals(points, faces)
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = zeros(Float64, N, N)
    D = zeros(Float64, N, N)
    H = zeros(Float64, N)

    for xkey in 1:N
        nx = normals[:, xkey]
        rx = points[:, xkey]

        for ykey in 1:N
            if xkey == ykey
                continue
            end

            ry = points[:, ykey]
            D[xkey, ykey] += dot(nx, ry-rx) / norm(ry-rx)^3 * deltaS[ykey]

        end
    end

    for xkey in 1:N
        H[xkey] = dot(H0, normals[:, xkey])
    end


    D_reg = D - Diagonal([sum(D[i, :]) for i in 1:N])
    return D_reg/(4pi) \ H
end


dHn = @time make_deltaH_normal(points, faces, 1., [0., 0., 1.])



fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")

(x, y, z) = [points[i,:] for i in 1:3]
(vx, vy, vz) = [(normals .* dHn')[i,:] . for i in 1:3]

ax[:scatter](x,y,z, s=2,color="k")
ax[:quiver](x,y,z,vx,vy,vz, length=30, arrow_length_ratio=0.5)

ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()
