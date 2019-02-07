using LinearAlgebra
using CSV
using JLD2
using ElTopo

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./stabilization.jl")
#using SurfaceGeometry

function normal_theor(vertex, coefs)
    normal = [2*vertex[1] / coefs.a^2, 2*vertex[2] / coefs.b^2, 2*vertex[3] / coefs.c^2]
    return normalize(normal)
end

function NormalsTheoretical(points, coef)
    all_theor_norm = zeros(Float64, 3, length(points[1,:]))
    for i in 1:size(points)[2]
        all_theor_norm[:, i] = normal_theor(points[:, i], coef)
    end
    return all_theor_norm
end


struct coefs
    a
    b
    c
end

#a, b = 1, 1
#c = 1/a/b
#coef = coefs(a, b, c)



function Normals(points, faces)
    all_normals = zeros(3, size(points,2))
    for i in 1:size(points)[2]
        mask = findall(x-> x == i, faces)
        neighbors = zeros(Int64, 3, length(mask))
        num_neighbors = length(mask)
        for n in 1:num_neighbors
            neighbors[:, n] = faces[:,mask[n][2]]
        end
        normal_i = zeros(3, 1)
        for j in 1:num_neighbors
            all_i = findall(neighbors .== i)
            v1 = neighbors[all_i[j][1], all_i[j][2]]
            v2 = neighbors[(all_i[j][1]) % 3 + 1, all_i[j][2]]
            v3 = neighbors[(all_i[j][1] + 1) % 3 + 1, all_i[j][2]]
            #println("v1 = $(norm(v1)), v2 = $(norm(v2)), v3 = $(norm(v3))")
            vec1 = points[:, v3] - points[:, v1]
            vec2 = points[:, v2] - points[:, v1]
            #println(norm(vec1), norm(vec2))
            normal_j = cross(vec1, vec2)
            angle_j = acos(dot(vec1, vec2) / norm(vec1) / norm(vec2))

            normal_i += angle_j * normal_j
        end
        normal_i = normalize(normal_i[:,1])
        all_normals[:, i] = normal_i
    end
    return all_normals
end

function NeighborVertices(vertex, faces)
    mask = findall(x-> x == vertex, faces)
    neighbors = zeros(Int64, 3, length(mask))
    num_neighbors = length(mask)
    for n in 1:num_neighbors
        neighbors[:, n] = faces[:,mask[n][2]]
    end
    return setdiff(neighbors, vertex)
end

function PotentialSimple(points,faces,hmag,H0;regularize=true,normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end


    A = zeros(Float64,size(points,2),size(points,2))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    for xkey in 1:size(points,2)

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


function HField(points,faces,psi)

    H = Array{Float64}(undef,size(points)...)
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        psix = psi[xkey]
        distances = Float64[]
        dphi = Float64[]

        for ykey in NeighborVertices(xkey,faces)
            y = points[:,ykey]
            distances = [distances; y-x]
            dphi = [dphi; psi[ykey]-psix]
        end
        A, B = distances, dphi
        A = transpose(reshape(A,3,div(length(A),3))) ### This looks unecesarry

        # linear least squares
        H[:,xkey] = inv(transpose(A)*A)*transpose(A)*B
    end
    return H
end


function HtField(points, faces, psi, normals)

    Ht = Array{Float64}(undef,3,size(points,2))
    H = HField(points,faces,psi)

    for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        P = I - nx*nx'
        Ht[:,xkey] = P*H[:,xkey]
    end

    return Ht
end

using FastGaussQuadrature
const NP = 10
const tt, ww = gausslegendre(NP)

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

function NormalFieldCurrent(points,faces,Ht,hmag,H0; eps=0.0001, normals=nothing)

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

    Hn = Array{Float64}(undef,size(points,2))

    for xkey in 1:size(points,2)

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

function InterfaceSpeedZinchenko(points,faces,forcen,etaP,gammap, normals)

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    velocityn = zeros(Float64,size(points,2))

    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]

        # vectorized
        fx = forcen[xkey]
        s = 0

        for ykey in 1:size(points,2)
            if ykey==xkey
                continue
            end

            y = points[:,ykey]
            ny = normals[:,ykey]

            # vectorized
            fy = forcen[ykey]

            ### I will need to check a missing 2
            s += vareas[ykey] * gammap ./8/pi/etaP * dot(y-x,nx+ny)/norm(y-x)^3 * (1-3*dot(y-x,nx)*dot(y-x,ny)/norm(y-x)^2)

            s += vareas[ykey]*1 ./8/pi/etaP * (dot(nx, ny)/norm(x-y) + dot(nx, y-x)*dot(ny, y-x)/norm(x-y)^3 ) * (fy - fx)
        end
        velocityn[xkey] = s
    end
    return velocityn
end




#
using Plots
# plots = Plo

# #
# pygui(true)
# pyplot()
#
# (x, y, z) = [zeros(1, size(points,2)) for i in 1:3]
# fig = figure(figsize=(7,7))
# ax = fig[:gca](projection="3d")
# for i in [50, 2000]
#     if i % 500 == 0
#     @load "./simul/points$i.jld2" points2
#     global x,y,z
#     (x, y, z) = [points2[j,:] for j in 1:3]
#     #end
#     ax[:scatter](x,y,z, s=3, label="$i")#,color="k")
# end
#ax[:quiver](x,y,z,vnx,vny,vnz, length=2, arrow_length_ratio=0.5)

#ax[:quiver](x,y,z,Hnx,Hny,Hnz, length=0.3, arrow_length_ratio=0.5)
#ax[:quiver](x,y,z,Htx,Hty,Htz, length=0.3, arrow_length_ratio=0.5)
#ax[:quiver](x,y,z,Hx,Hy,Hz, length=0.3)
#
# ax[:set_xlim](-2,2)
# ax[:set_ylim](-2,2)
# ax[:set_zlim](-2,2)
# ax[:set_xlabel]("x axis")
# ax[:set_ylabel]("y axis")
# ax[:set_zlabel]("z axis")
# ax[:legend]()
# fig[:show]()
