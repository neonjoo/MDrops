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


function PotentialSimple_par(points,faces,normals,mu,H0;regularize=true)
    # return the magnetic potential psi on the surface of teh drop

    hmag = mu
    A = zeros(Float64,size(points,2),size(points,2))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    Threads.@threads for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]

        for ykey in 1:size(points,2)
            if xkey==ykey
                continue
            end
            y = points[:,ykey]
            ny = normals[:,ykey]

            A[ykey,xkey] = dot(y-x,ny)/norm(y-x)^3 * vareas[ykey]
        end
    end

    B = zeros(Float64,size(points,2))
    for xkey in 1:size(points,2)
        B[xkey] = 2*dot(H0,points[:,xkey])/(hmag+1)
    end

    if regularize==true
        A = A'
        reg_A = A - Diagonal(Float64[sum(A[i,:]) for i in 1:size(A,2)])
        psi = (I * (1- (hmag-1)/(hmag+1)) - 1/2/pi * (hmag-1)/(hmag+1) * reg_A) \ B
    else
        A = A*(hmag-1)/(hmag+1)/2/pi
        A = A'
        psi = (I - A)\B
    end

    return psi
end

function HtField_par(points, faces, psi, normals)
    # return the tangential component of the magnetic field on the droplet surface
    # returns a vector
    Ht = Array{Float64}(undef,3,size(points,2))
    H = HField(points,faces,psi)

    Threads.@threads for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        P = I - nx*nx'
        Ht[:,xkey] = P*H[:,xkey]
    end

    return Ht
end

function NormalFieldCurrent_par(points,faces,normals,Ht,mu,H0; eps=0.0001)
    # return normal component of magnetic field on the surface of the droplet
    # returns a scalar

    hmag = mu

    for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        #Ht[:,xkey] = (eye(3) - nx*nx')*Ht[:,xkey]
    end

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end


    function qs(xi,eta,v1,v2,v3,x,nx,Htx)
        y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
        Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
        ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
        s = - dot(nx,cross(Hty - Htx,cross(ny,-(y-x)/norm(y-x)^2)))
    end

    NP = 10
    tt, ww = gausslegendre(NP)

    function strquad(q::Function,x1,x2,x3)

        B = dot(x3-x1,x2-x1)/norm(x2-x1)^2
        C = norm(x3-x1)^2/norm(x2-x1)^2
        hS = norm(cross(x2-x1,x3-x1))

        s = 0
        for i in 1:NP
            Chi = pi/4*(1 + tt[i])

            R = 1/(cos(Chi) + sin(Chi))
            si = 0
            for j in 1:NP
                rho = R/2*(1 + tt[j])
                si += q(rho*cos(Chi),rho*sin(Chi))*ww[j]
            end

            s += si*R/2 / sqrt(cos(Chi)^2 + B*sin(2*Chi) + C*sin(Chi)^2) * ww[i]
        end
        s *= pi/4

        return s*hS/norm(x2 - x1)
    end

    Hn = Array{Float64}(undef,size(points,2))

    Threads.@threads for xkey in 1:size(points,2)

        nx = normals[:,xkey]
        x = points[:,xkey] + eps*nx
        Htx = Ht[:,xkey]

        s = 0
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            y = points[:,ykey]
            ny = normals[:,ykey]
            Hty = Ht[:,ykey]

            #s += dot(nx,-(Hty-Htx)*dot((y-x)/norm(y-x)^3,ny)) * vareas[ykey]
            s += dot(nx,-(Hty-Htx)*dot((y-x)/norm(y-x)^3,ny)) * vareas[ykey]
            s += -dot(nx,cross(Hty-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * vareas[ykey]
            #s += dot(nx,cross(cross(ny,),-(y-x)/norm(y-x)^3)) * vareas[ykey]
        end

        ### Making a proper hole
        for (v2,v3) in DoubleVertexVRing(xkey,faces)
            area = norm(cross(points[:,v2]-x,points[:,v3]-x))/2

            ny = normals[:,v2]
            y = points[:,v2]
            s -= -dot(nx,cross(Ht[:,v2]-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * area/3

            ny = normals[:,v3]
            y = points[:,v3]
            s -= -dot(nx,cross(Ht[:,v3]-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * area/3

            ### Singular triangle integration

            #s += strquad((xi,eta) -> qs(xi,eta,xkey,v2,v3,x,nx,Htx),x,points[:,v2],points[:,v3],abstol=abs(s/100))[1]
            s += strquad((xi,eta) -> qs(xi,eta,xkey,v2,v3,x,nx,Htx),x,points[:,v2],points[:,v3])
        end

        #println("xkey is $xkey")

        Hn[xkey] = dot(H0,nx)/hmag + 1/4/pi * (1-hmag)/hmag * s
    end

    return Hn
end

function make_magvelocities_par(vertices, normals, lambda, Bm, mu, Hn_2, Ht_2)
    #without the use of levi civita tensor
    # lambda = int viscosity / ext viscosity
    # Returns vertex velocities from the Stokes flow integral equations
    # using Wielandt's deflation as describeed in Zinchenko 1997, Das 2017
    # w + T*w + 4*pi*C*w + B*w = F
    println("calculating velocities")
    deltaS = make_dS(points,faces)

    #fbar = mu*(mu-1)/8pi * Hfieldn.^2 + (mu-1)/8pi * Hfieldt.^2
    fbar = mu*(mu-1)/8pi * Hn_2 + (mu-1)/8pi * Ht_2


    S = sum(deltaS)
    xc = [0.,0.,0.]
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        xc = xc + rx * deltaS[x] / S
    end

    F = zeros(3*size(vertices,2),1)
    T = zeros(3*size(vertices,2),3*size(vertices,2))
    #C = C1+C2
    C1 = zeros(3*size(vertices,2),3*size(vertices,2))
    C2 = zeros(3*size(vertices,2),3*size(vertices,2))
    B = zeros(3*size(vertices,2),3*size(vertices,2))

    M = zeros(3,3)
    # fill M
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        # Matrix(1.0I, 3, 3) <- 3x3 Identity matrix
        M += (Matrix(1.0I, 3, 3) * dot(rx - xc, rx - xc) -
            (rx - xc) * (rx - xc)') * deltaS[x]
            #   ^ changed the transpose
    end
    Minv = inv(M) # calculate inverse of M

    Threads.@threads for y = 1:size(vertices,2)
        ry = vertices[:,y]

        for j = 1:3

            for x = 1:size(vertices,2)
                rx = vertices[:,x]

                Gji = (Matrix(1.0I, (3,3)) / norm(rx-ry) +
                (rx-ry) * (rx-ry)' / norm(rx-ry)^3)

                if x != y
                   F[3*y-3+j] +=
                       1/(8*pi) *(
                       dot(rx-ry,normals[:,x]) * normals[j,y]+
                       dot(rx-ry,normals[:,y]) * normals[j,x]+
                       (1-dot(normals[:,x],normals[:,y])) * (rx[j]-ry[j])-
                       3*(rx[j]-ry[j])/norm(rx-ry)^2*
                       dot(normals[:,x]+normals[:,y],rx-ry)*
                       dot(rx-ry,normals[:,x])) * deltaS[x]/norm(rx-ry)^3

                # magnetic part
                   F[3*y-3+j] += Bm/(8*pi) *
                   (fbar[x] - fbar[y]) *
                   dot(Gji[j,:], normals[:,x]) * deltaS[x]
                #else
                #polar integration scheme
                end

                for i = 1:3
                    if x != y
                       T[3*y-3+j,3*x-3+i] = -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]

                       T[3*y-3+j,3*y-3+i] -= -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]
                    end

                    B[3*y-3+j,3*x-3+i] = -4*pi/S * normals[j,y] * normals[i,x] * deltaS[x]

                    C1[3*y-3+j,3*x-3+i] =  deltaS[x]/S

                    for k = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            Minv[k,j]*(rx[k] - xc[k])*deltaS[x]*(ry[i] - xc[i]) -
                            Minv[i,j]*(rx[k] - xc[k])*deltaS[x]*(ry[k] - xc[k]) -
                            Minv[k,k]*(rx[j] - xc[j])*deltaS[x]*(ry[i] - xc[i]) +
                            Minv[i,k]*(rx[j] - xc[j])*deltaS[x]*(ry[k] - xc[k])
                    end
                    if i==j
                        for k = 1:3
                        for m = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            Minv[k,k]*(rx[m] - xc[m])*deltaS[x]*(ry[m] - xc[m]) -
                            Minv[k,m]*(rx[k] - xc[k])*deltaS[x]*(ry[m] - xc[m])
                        end
                        end
                    end

                end
            end

        end
    end
    C = C1 + C2

    # w + T*w + 4*pi*C*w + B*w = F
    #w = (eye(size(T)) + T + 4*pi*C + B) \ F; old and bad
    w = (Matrix(1.0I, size(T)) + T + (lambda-1)/2*C + (lambda-1)/(8*pi)*B) \ F

    # w' = C*w
    # v = w + (lambda-1)/2 * w'
    varr = w + (lambda-1)/2 * C*w


    velocities = reshape(varr ,3,size(vertices,2))
    return velocities
end


# @load "./meshes/points_critical_hyst_2_21.jld2"
# @load "./meshes/faces_critical_hyst_2_21.jld2"
#
# points = Array{Float64}(points')
# faces = Array{Int64}(faces')

# a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])
#

# edges = make_edges(faces)
# connectivity = make_connectivity(edges)
# normals = Normals(points, faces)
# (normals, CDE) = make_normals_spline(points, connectivity, edges, normals)

# # #
# # ## setting simulation parameters
# H0 = [0., 0., 1.]
# mu = 10.
# lambda = 10.
# Bm = 3.

## calculating the magnetic field on the surface
# @time psi = PotentialSimple(points, faces, normals, mu, H0)
# @time psi_par = PotentialSimple_par(points, faces, normals, mu, H0)
#
# Ht_vec = HtField_par(points, faces, psi, normals) # a vector
# # Ht = sqrt.(sum(Ht_vec.^2,dims=1))'
# @time Hn_norms = NormalFieldCurrent_par(points, faces, normals, Ht_vec, mu, H0) # a scalar
# Hn = normals .* Hn_norms'
#
# Hn_2 = sum(Hn.^2, dims=1)
# Ht_2 = sum(Ht_vec.^2, dims=1)
#
# @time v = make_magvelocities(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
# @time vp = make_magvelocities_par(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
#
# vvecs = v
#
# @time v_stab = make_Vvecs_conjgrad(normals, faces, points, vvecs, 1e-6, 500)
# @time v_stabp = make_Vvecs_conjgrad_par(normals, faces, points, vvecs, 1e-6, 500)
