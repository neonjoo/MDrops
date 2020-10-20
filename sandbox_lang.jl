using Pkg
pkg"activate ."
pkg"resolve"

using JLD2
using StatsBase
using LinearAlgebra
using FastGaussQuadrature
using Optim

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")

## making the mesh
points, faces = expand_icosamesh(;R=1,depth=3)
points = Array{Float64}(points)
faces = Array{Int64}(faces)




function make_F_par2(triangles, vertices, V)

    Ntriangles = size(triangles, 2)
    F_values = zeros(Ntriangles)

    Threads.@threads for i = 1:Ntriangles
        x1 = vertices[:,triangles[1,i]]
        x2 = vertices[:,triangles[2,i]]
        x3 = vertices[:,triangles[3,i]]

        v1 = V[:,triangles[1,i]]
        v2 = V[:,triangles[2,i]]
        v3 = V[:,triangles[3,i]]

        a = norm(x2 - x1)
        b = norm(x3 - x2)
        c = norm(x1 - x3)

        Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

        A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
            (4*Cdelta * (a^2 + b^2 + c^2)^3)
        B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
            (4*Cdelta * (a^2 + b^2 + c^2)^3)
        C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
            (4*Cdelta * (a^2 + b^2 + c^2)^3);

        dCdeltadt = -A * dot(x2 - x1, v2 - v1) +
                    -B * dot(x3 - x2, v3 - v2) +
                    -C * dot(x1 - x3, v1 - v3)

        F_values[i] = 0.4 / Cdelta^2 * dCdeltadt^2 +
            2*( dot(x2 - x1, v2 - v1))^2 +
            2*( dot(x3 - x2, v3 - v2))^2 +
            2*( dot(x1 - x3, v1 - v3))^2
    end

    return sum(F_values)
end


function make_tanggradF_par2(normals,triangles, vertices, V,neighbor_triang_masks)

        Nvertices = size(vertices, 2)
        gradF = zeros(3,Nvertices)

        Threads.@threads for i = 1:Nvertices

            # finds the i-th triangle indices in the triangle 2D array
            mask = neighbor_triang_masks[i]
            num = length(mask)
            (row, col) = zeros(Int64, 1, num), zeros(Int64, 1, num)
            for n in 1:num
                row[n], col[n] = mask[n][1], mask[n][2]
            end


            for this_triang = 1:length(row)
                j1 = (row[this_triang])%(3) + 1 # from 1 to 3
                j2 = (row[this_triang] + 1)%(3) +1 # from 1 to 3

                x1 = vertices[:,triangles[row[this_triang],col[this_triang]]]
                x2 = vertices[:,triangles[j1, col[this_triang]]]
                x3 = vertices[:,triangles[j2, col[this_triang]]]

                v1 = V[:, triangles[row[this_triang],col[this_triang]]]
                v2 = V[:, triangles[j1, col[this_triang]]]
                v3 = V[:, triangles[j2, col[this_triang]]]

                a = norm(x2 - x1)
                b = norm(x3 - x2)
                c = norm(x1 - x3)

                Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

                A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)
                B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)
                C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)

                dCdeltadt = -A * dot(x2 - x1, v2 - v1) +
                            -B * dot(x3 - x2, v3 - v2) +
                            -C * dot(x1 - x3, v1 - v3)

                t1 = 0.4 / Cdelta^2 * 2 * dCdeltadt * ( A*(x2 - x1) + C*(x3 - x1))
                t2 = -4*dot(x2 - x1, v2 - v1) * (x2 - x1)
                t3 = -4*dot(x3 - x1, v3 - v1) * (x3 - x1)

                gradF[:,i] = gradF[:,i] + t1 + t2 + t3

            end

            tang_proj = I - normals[:,i] * normals[:,i]'
            gradF[:,i] = tang_proj * gradF[:,i]

        end

        return gradF

end

function make_gradF_par2(normals,triangles, vertices, V,neighbor_triang_masks)

        Nvertices = size(vertices, 2)
        gradF = zeros(3,Nvertices)

        Threads.@threads for i = 1:Nvertices

            # finds the i-th triangle indices in the triangle 2D array
            mask = neighbor_triang_masks[i]
            num = length(mask)
            (row, col) = zeros(Int64, 1, num), zeros(Int64, 1, num)
            for n in 1:num
                row[n], col[n] = mask[n][1], mask[n][2]
            end


            for this_triang = 1:length(row)
                j1 = (row[this_triang])%(3) + 1 # from 1 to 3
                j2 = (row[this_triang] + 1)%(3) +1 # from 1 to 3

                x1 = vertices[:,triangles[row[this_triang],col[this_triang]]]
                x2 = vertices[:,triangles[j1, col[this_triang]]]
                x3 = vertices[:,triangles[j2, col[this_triang]]]

                v1 = V[:, triangles[row[this_triang],col[this_triang]]]
                v2 = V[:, triangles[j1, col[this_triang]]]
                v3 = V[:, triangles[j2, col[this_triang]]]

                a = norm(x2 - x1)
                b = norm(x3 - x2)
                c = norm(x1 - x3)

                Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

                A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)
                B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)
                C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)

                dCdeltadt = -A * dot(x2 - x1, v2 - v1) +
                            -B * dot(x3 - x2, v3 - v2) +
                            -C * dot(x1 - x3, v1 - v3)

                t1 = 0.4 / Cdelta^2 * 2 * dCdeltadt * ( A*(x2 - x1) + C*(x3 - x1))
                t2 = -4*dot(x2 - x1, v2 - v1) * (x2 - x1)
                t3 = -4*dot(x3 - x1, v3 - v1) * (x3 - x1)

                gradF[:,i] = gradF[:,i] + t1 + t2 + t3

            end

        end

        return gradF

end

function make_Vvecs_conjgrad_par2(normals,triangles, vertices, vvecs, epsilon, maxIters)

    println("passive stabbing")
    # finds indices of triangles near point i
    neighbor_triang_masks = [findall(x-> x == i, triangles) for i in 1:size(vertices, 2)]
    # first gradient descent
    f = make_tanggradF_par2(normals,triangles, vertices, vvecs, neighbor_triang_masks)
    gradFv = make_gradF_par2(normals, triangles, vertices, vvecs, neighbor_triang_masks)
    gradFf = make_gradF_par2(normals, triangles, vertices, f, neighbor_triang_masks)

    ksi = - sum(sum(gradFv .* f, dims=2)) / sum(sum(gradFf .* f, dims=2))

    V = vvecs + ksi*f
    Vp = vvecs

    F = make_F_par2(triangles, vertices, V)

    # then conjugated gradient()
    for i = 1:maxIters
        f = make_tanggradF_par2(normals,triangles, vertices, V, neighbor_triang_masks)
        gradFv = make_gradF_par2(normals, triangles, vertices, V, neighbor_triang_masks)
        gradFvp = make_gradF_par2(normals, triangles, vertices, Vp, neighbor_triang_masks)
        gradFf = make_gradF_par2(normals, triangles, vertices, f, neighbor_triang_masks)
        Ff = make_F_par2(triangles, vertices, f)
        Fdeltav = make_F_par2(triangles, vertices, V-Vp)

        a1 = sum(sum(gradFv .* f, dims=2))
        b1 = Ff
        c1 = sum(sum((gradFv .- gradFvp) .* f, dims=2))
        a2 = sum(sum(gradFv .* (V.-Vp), dims=2))
        b2 = sum(sum(gradFf .* (V.-Vp), dims=2))
        c2 = Fdeltav

        ksi = (a2*c1 - 2*a1*c2) / (4*b1*c2 - b2*c1)
        eta = (2*a2*b1 - a1*b2) / (b2*c1 - 4*b1*c2)

        Vtemp = V
        V = V .+ ksi*f .+ eta*(V.-Vp)
        Vp = Vtemp

        Fp = F
        F = make_F_par2(triangles, vertices, V)

        if (Fp - F)/F < epsilon
            println("improved tangential velocities")
            break
        end
        if i == maxIters
            println("tangential velocities not fully converged")
        end
    end

    return V
end




edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
#(normals, CDE, AB) = make_normals_parab(points, connectivity, edges, normals)

# #
# ## setting simulation parameters
H0 = [0., 0., 1.]
mu = 10.
lambda = 10.
Bm = 3.

# calculating the magnetic field on the surface
@time psi = PotentialSimple(points, faces, normals, mu, H0)
@time psi_par = PotentialSimple_par(points, faces, normals, mu, H0)

Ht_vec = HtField_par(points, faces, psi, normals) # a vector
# Ht = sqrt.(sum(Ht_vec.^2,dims=1))'
@time Hn_norms = NormalFieldCurrent_par(points, faces, normals, Ht_vec, mu, H0) # a scalar
Hn = normals .* Hn_norms'

Hn_2 = sum(Hn.^2, dims=1)
Ht_2 = sum(Ht_vec.^2, dims=1)

@time v_phys = make_magvelocities_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)

vvecs = v_phys

println(333333)
@time v_stab = make_Vvecs_conjgrad(normals, faces, points, vvecs, 1e-6, 500)
@time v_stab_par = make_Vvecs_conjgrad_par(normals, faces, points, vvecs, 1e-6, 500)
@time v_stab_par2 = make_Vvecs_conjgrad_par2(normals, faces, points, vvecs, 1e-6, 500)
