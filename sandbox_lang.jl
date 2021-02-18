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

function make_gradE_par(points,faces,closefaces,hsq; p=50.,r=100.)
    Cet = 0.144337567297 # sqrt(3) / 12 # target Cdelta value
    Nvertices = size(points, 2)
    gradE = zeros(3,Nvertices)
    for i = 1:Nvertices
        for faceind in closefaces[:,i]
            if faceind == 0
                break
            end
        # order radius vectors x so that first one points to vertex i
        iind = findfirst(x->x==i,faces[:,faceind])
        x1 = points[:,faces[iind,faceind]]
        x2 = points[:,faces[iind%3+1,faceind]]
        x3 = points[:,faces[(iind+1)%3+1,faceind]]

        hsq1 = hsq[faces[iind,faceind]]
        hsq2 = hsq[faces[iind%3+1,faceind]]
        hsq3 = hsq[faces[(iind+1)%3+1,faceind]]

        a = norm(x2 - x1)
        b = norm(x3 - x2)
        c = norm(x1 - x3)

        Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        x12 = x2 - x1
        x13 = x3 - x1
        h12sq = 0.5*(hsq1 + hsq2)
        h13sq = 0.5*(hsq1 + hsq3)

        dCdx = 1/(4*Cdelta*(a^2 + b^2 + c^2)^3) * (
                -(a^4 + b^4 + c^4 - (a^2 + b^2 + c^2)*a^2) * x12
                -(a^4 + b^4 + c^4 - (a^2 + b^2 + c^2)*c^2) * x13
                )

        gradE[:,i] = gradE[:,i] +
            r*(Cet/Cdelta)^(r-1) * Cet * (-1/Cdelta^2) * dCdx +
            0.5*p*( 0.5*( dot(x12,x12)/h12sq + h12sq/dot(x12,x12)) )^(p-1) *(-1)* (1/h12sq - h12sq/dot(x12,x12)^2) * x12 +
            0.5*p*( 0.5*( dot(x13,x13)/h13sq + h13sq/dot(x13,x13)) )^(p-1) *(-1)* (1/h13sq - h13sq/dot(x13,x13)^2) * x13
        end
    end
    return gradE
end


function active_stabilize_old_surface_par(points_old,CDE_old,normals_old,points0::Array{Float64,2},faces::Array{Int64,2},connectivity::Array{Int64,2},edges::Array{Int64,2};
    deltakoef=0.01, R0=1.0, gamma=0.25, p=50., r=100., checkiters=100, maxiters=1000,critSc = 0.75,critCdelta = 1.15)
    # actively rearange vertices on a surfaces given by fitted paraboloids
    # this has been modified to use points before splitting for surface determination
    # as per Zinchenko(2013)
    println("active stabilization")
    points = copy(points0)

    closefaces = make_closefaces(faces)
    dS = make_dS(points,faces)

    k1 = Array{Float64}(undef, size(points,2))
    k2 = Array{Float64}(undef, size(points,2))
    for i = 1:size(points,2)
        r0 = points[:,i]
        minind = argmin(sum((points_old .- r0).^2,dims=1))[2]
        r0 = to_local(r0-points_old[:,minind],normals_old[:,minind])
        k1[i],k2[i] =  make_pc_local(CDE_old[:,minind],r0[1],r0[2])
    end

    LAMBDA = k1.^2 + k2.^2 .+ 0.004/R0^2
    K = 4/(sqrt(3) * size(faces,2)) * sum(LAMBDA.^gamma .* dS)
    hsq = K * LAMBDA.^(-gamma)

    no_improvement = true
    # initiallize improvement criteria
    xij = make_edge_lens(points,edges)
    hij = sqrt.(0.5*(hsq[edges[1,:]].^2 + hsq[edges[2,:]].^2))

    Sc0 = maximum(xij./hij) / minimum(xij./hij)
    Cdelta_min0 = minimum(make_Cdeltas(points, faces))
    for iter = 1:maxiters
        #println(iter)
        gradE = make_gradE_par(points,faces,closefaces,hsq; p=p,r=r)
        delta = deltakoef * minimum(make_min_edges(points,connectivity) ./ sum(sqrt.(gradE.^2),dims=1))

        #println("dPoints= ",delta*maximum(sum(sqrt.(gradE.^2),dims=1)))
        #println("E = ", make_E(points,faces,hsq; p=p,r=r))


        points = points - delta * gradE

        #project points on the drop
        Threads.@threads for i = 1:size(points,2)
            points[:,i] = project_on_drop(points_old,CDE_old,normals_old,points[:,i])
        end

        # recalculate the mesh parameters
        dS = make_dS(points,faces)
        #k1,k2 = make_pc(CDE)
        Threads.@threads for i = 1:size(points,2)
            r0 = points[:,i]
            minind = argmin(sum((points_old .- r0).^2,dims=1))[2]
            r0 = to_local(r0-points_old[:,minind],normals_old[:,minind])
            k1[i],k2[i] =  make_pc_local(CDE_old[:,minind],r0[1],r0[2])
        end
        LAMBDA = k1.^2 + k2.^2 .+ 0.004/R0^2
        K = 4/(sqrt(3) * size(faces,2)) * sum(LAMBDA.^gamma .* dS)
        hsq = K * LAMBDA.^(-gamma)
        if iter < checkiters
            xij = make_edge_lens(points,edges)
            hij = sqrt.(0.5*(hsq[edges[1,:]].^2 + hsq[edges[2,:]].^2))

            Sc = maximum(xij./hij) / minimum(xij./hij)
            Cdelta_min = minimum(make_Cdeltas(points, faces))

            #println("Sc/Sc0 = ",Sc/Sc0)
            #println("Cdelta/Cdelta0 = ",Cdelta_min/Cdelta_min0)
            if Sc > critSc*Sc0 && Cdelta_min < critCdelta*Cdelta_min0
                no_improvement = false
            end
        end

        if iter == checkiters
            if no_improvement == true
                println("no significant improvement achieved")
                println("reversing changes")
                points = points0
                break
            else
                println("improvement detected in the first ", checkiters, " iterations")
                println("iterating for ", maxiters - checkiters, " more iterations")
            end
        end

        # e = make_E(points, faces, hsq,p=p, r=r)
        # println("E = $e")

        if iter%500 == 0
            println("iteration ",iter)

        end
    end
    return points
end



## making the mesh
points, faces = expand_icosamesh(;R=1,depth=2)
points = Array{Float64}(points)
faces = Array{Int64}(faces)




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
@time psi_par = PotentialSimple_par(points, faces, normals, mu, H0)

Ht_vec = HtField_par(points, faces, psi, normals) # a vector
# Ht = sqrt.(sum(Ht_vec.^2,dims=1))'
@time Hn_norms = NormalFieldCurrent_par(points, faces, normals, Ht_vec, mu, H0) # a scalar
Hn = normals .* Hn_norms'
normals, CDE, AB = make_normals_parab(points, connectivity, normals)
Hn_2 = sum(Hn.^2, dims=1)
Ht_2 = sum(Ht_vec.^2, dims=1)

@time v_phys = make_magvelocities_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)

vvecs = v_phys
points_new = points + 0.05 * v_phys

normals_new, CDE_new, AB = make_normals_parab(points, connectivity, normals)


@time Juno.@profiler p_stab = active_stabilize_old_surface(points,
                        CDE, normals, points_new,
                        faces, connectivity, edges)

@time Juno.@profiler p_stab_par = active_stabilize_old_surface_par(points,
                        CDE, normals, points_new,
                        faces, connectivity, edges)

# println(333333)
# @time v_stab = make_Vvecs_conjgrad(normals, faces, points, vvecs, 1e-6, 500)
# @time v_stab_par = make_Vvecs_conjgrad_par(normals, faces, points, vvecs, 1e-6, 500)
# @time v_stab_par2 = make_Vvecs_conjgrad_par2(normals, faces, points, vvecs, 1e-6, 500)
