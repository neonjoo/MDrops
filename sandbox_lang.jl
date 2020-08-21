using Pkg
pkg"activate ."
pkg"resolve"

using Plots
using StatsBase

#p = Plots
#p.pyplot()
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")

function gauss_nonsingular(f::Function, r1,r2,r3,gaussorder)
    #f - function of r
    u, wu = gausslegendre(gaussorder) # u from -1 to 1
    v, wv = u, wu # for convenient notation

    e_xi = r2-r1
    e_eta = r3-r1
    hs = norm(cross(e_xi,e_eta))

    function make_r(u,v,r1,r2,r3)
        xi = (1+u)/2
        eta = (1-u)*(1+v)/4
        zeta = 1-xi-eta

        r = zeta*r1 + xi*r2 + eta*r3

        return r
    end

    for i = 1:length(u)
        for k = 1:length(v)
            if k == 1 && i==1
                # initializing here to avoid guessing the type of output of f()
                intval = wu[i]*(1-u[i]) * wv[k]*f( make_r(u[i],v[k],r1,r2,r3) )
            else
            intval += wu[i]*(1-u[i]) * wv[k]*f( make_r(u[i],v[k],r1,r2,r3) )
            end
        end
    end
    intval *= hs/8

    return intval
end

function gauss_weaksingular(q::Function,r1,r2,r3,gaussorder)
    # integral of q(r)/|r-r1|, where q(r) is nonsingular
    # singularity on r1

    u, wu = gausslegendre(gaussorder) # u from -1 to 1
    v, wv = u, wu # for convenient notation

    B = dot(r3-r1,r2-r1)/norm(r2-r1)^2
    C = norm(r3-r1)^2/norm(r2-r1)^2

    e_xi = r2-r1
    e_eta = r3-r1
    hs = norm(cross(e_xi,e_eta))

    function make_R(u)
        chi = (u+1)/4 * pi
        return 1/(cos(chi) + sin(chi))
    end

    function make_r(u,v,r1,r2,r3)
        chi = (u+1)/4 * pi
        rho = (v+1)/2 * make_R(u)

        xi = rho*cos(chi)
        eta = rho*sin(chi)

        zeta = 1-xi-eta
        r = zeta*r1 + xi*r2 + eta*r3
        return r
    end

    function make_chi(u)
        return (u+1)/4 * pi
    end

    for i = 1:length(u)
        chi = make_chi(u[i])
        for k = 1:length(v)
            if k == 1 && i==1
                # initializing here to avoid guessing the type of output of f()
                intval = wu[i]*make_R(u[i])/sqrt( cos(chi)^2 + B*sin(2*chi) + C*sin(chi)^2 ) *
                 wv[k]*q( make_r(u[i],v[k],r1,r2,r3) )
            else
                intval += wu[i]*make_R(u[i])/sqrt( cos(chi)^2 + B*sin(2*chi) + C*sin(chi)^2 ) *
                 wv[k]*q( make_r(u[i],v[k],r1,r2,r3) )
            end
        end
    end
    intval *= pi*hs / (8 * norm(r2-r1))

    return intval

end

points, faces = expand_icosamesh(R=1, depth=2)
points = Array{Float64}(points)
faces = Array{Int64}(faces)
normals = Normals(points, faces)
H0 = [0., 0., 1.]

function make_L_curved(points, faces)
    normals = Normals(points, faces)
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = zeros(Float64, N)

    edges = make_edges(faces)
    connectivity = make_connectivity(edges)
    normals, CDE = make_normals_spline(points, connectivity, edges, normals)
    curvas = make_pc(CDE)
    divN = curvas[1] + curvas[2]
    for ykey in 1:N
        ny = normals[:, ykey]
        ry = points[:, ykey]

        for xkey in 1:N

            if ykey == xkey
                continue
            end
            nx = normals[:, xkey]
            rx = points[:, xkey]

            int = dot(ry-rx, nx)*(nx-ny) / norm(ry-rx)^3 + nx * divN[xkey] / norm(ry-rx)

            L[ykey] += dot(ny, int) * deltaS[xkey] / 4pi
        end
    end

    return L
end

function make_L_curved_gauss(points, faces; gaussorder=3)
    normals = Normals(points, faces)
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = zeros(Float64, N)

    edges = make_edges(faces)
    connectivity = make_connectivity(edges)
    normals, CDE = make_normals_spline(points, connectivity, edges, normals)
    curvas = make_pc(CDE)
    divN = curvas[1] + curvas[2]

    for ykey in 1:N
        ny = normals[:, ykey]
        ry = points[:, ykey]

        function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
            A = [x1 x2 x3] # matrix of vertex radiusvecotrs
            B = [n1 n2 n3] # matrix of vertex normals

            zeta_xi_eta = A \ x # find local triangle parameters

            n = B * zeta_xi_eta
            return n/norm(n)
        end

        function make_divn(x,x1,x2,x3,divn1,divn2,divn3) # normal linear interpolation
            A = [x1 x2 x3] # matrix of vertex radiusvecotrs
            B = [divn1 divn2 divn3] # matrix of vertex curvatures

            zeta_xi_eta = A \ x # find local triangle parameters

            divn = B * zeta_xi_eta
            return divn[1]
        end

        function make_L(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3; y=points[:,ykey],ny=normals[:,ykey])
            nx = make_n(x,x1,x2,x3,n1,n2,n3)
            divnx = make_divn(x, x1, x2, x3, divn1, divn2, divn3)
            r = y - x

            int = dot(r, nx)*(nx-ny) / norm(r)^3 + nx * divnx/ norm(r)
            return dot(ny, int) / 4pi
        end

        for i in 1:size(faces, 2)
            x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
            n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
            divn1, divn2, divn3 = [divN[faces[j,i]] for j in 1:3]

            #int = dot(ry-rx, nx)*(nx-ny) / norm(ry-rx)^3 + nx * divN[xkey] / norm(ry-rx)
            #L[ykey] += dot(ny, int) * deltaS[xkey] / 4pi

            L[ykey] += gauss_nonsingular(x->make_L(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3), x1,x2,x3,gaussorder)
        end
    end

    return L
end

function make_L_sing(points, faces; gaussorder=5)
    normals = Normals(points, faces)
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = zeros(Float64, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end


    for ykey in 1:N
        function make_L(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
            nx = make_n(x,x1,x2,x3,n1,n2,n3)
            r = y - x

            return dot(r, nx-ny) / norm(r)^3 / 4pi
        end

        function make_L_times_norm_r(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey], ny=normals[:,ykey])
            nx = make_n(x,x1,x2,x3,n1,n2,n3)
            r = y - x

            return dot(r, nx-ny) / norm(r)^2 / 4pi
        end

        ny = normals[:, ykey]
        ry = points[:, ykey]

        for i in 1:size(faces, 2)
            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

                L[ykey] += gauss_nonsingular(x->make_L(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 1) % 3 + 1,i]]
                #println("L = $L,eeee $(L[ykey])")
                #L[ykey] +=
                L[ykey] += gauss_weaksingular(x->make_L_times_norm_r(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #println("$Ltemp, $(L[ykey]), $ykey")
            end

        end
    end

    return -L # implement the minus in the end. For a sphere L=-1 for all points
end

function make_deltaH_normal(points, faces, alpha, H0)
    normals = Normals(points, faces)
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = make_L_sing(points, faces; gaussorder=5)
    L = Diagonal(alpha/(alpha-1) .- L) # alpha term should be positive
    F = zeros(Float64, N, N) # integral equation matrix
    H = zeros(Float64, N)

    for ykey in 1:N
        ny = normals[:, ykey]
        ry = points[:, ykey]

        for xkey in 1:N
            if xkey == ykey
                continue
            end
            nx = normals[:, xkey]
            rx = points[:, xkey]
            F[ykey, xkey] +=  dot(ny, ry-rx) / norm(ry-rx)^3 * deltaS[xkey] / (4pi)

        end
    end

    for ykey in 1:N
        H[ykey] = dot(H0, normals[:, ykey]) # no minus
    end

    F_reg = F - Diagonal([sum(F[i, :]) for i in 1:N])
    deltaH_normal = (L - F_reg) \ H

    return deltaH_normal' / (alpha - 1)
end

#points_csv= CSV.read("./meshes/points_ellipse_fewN.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_ellipse_fewN.csv", header=0)
# 
# points = convert(Array, points_csv)
# faces = convert(Array, faces_csv)
#
# dHn = @time make_deltaH_normal(points, faces, 2, [0., 0., 1.])
#
# Hs2_e_L5 = make_deltaH_normal(points, faces, 2, [0., 0., 1.])
# Hs5_el = make_deltaH_normal(points, faces, 5, [0., 0., 1.])
# Hs30 = make_deltaH_normal(points, faces, 30, [0., 0., 1.])
#
# L_kurva3 = @time make_L_curved_gauss(points, faces; gaussorder=3)
# L_kurva5 = @time make_L_curved_gauss(points, faces; gaussorder=5)
# #L_kurva10 = @time make_L_curved_gauss(points, faces; gaussorder=10)
# L_sing3 = @time make_L_sing(points, faces; gaussorder=3)
# L_sing5 = @time make_L_sing(points, faces; gaussorder=5)
# L_sing10 = @time make_L_sing(points, faces; gaussorder=10)
#
# normals = Normals(points, faces)
# edges = make_edges(faces)
# connectivity = make_connectivity(edges)
# #normals, CDE = make_normals_spline(points, connectivity, edges, normals)
#
# mu=2
# psi = PotentialSimple(points, faces, mu, H0; normals = normals)
# Ht = HtField(points, faces, psi, normals)
# Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
# Hn = normals .* Hn_norms'


# using PyPlot

# fig = figure(figsize=(7,7))
# ax = fig[:gca](projection="3d")
# (x, y, z) = [points[i,:] for i in 1:3]
# (vx, vy, vz) = [(normals .* Ht2)[i,:] for i in 1:3]
# #ax[:quiver](x,y,z,vx,vy,vz, length=0.6, arrow_length_ratio=0.5, color=:black)
# (vx, vy, vz) = [(normals .* Hs2)[i,:] for i in 1:3]
# #ax[:quiver](x,y,z,vx,vy,vz, length=0.6, arrow_length_ratio=0.5, color=:blue)
# (vx, vy, vz) = [(normals .* erd_norm)[i,:] for i in 1:3]
# ax[:quiver](x,y,z,vx,vy,vz, length=0.6, arrow_length_ratio=0.5, color=:red)
# ax[:scatter](x,y,z, s=2,color="k")
# #ax[:quiver](x,y,z,vx,vy,vz, length=50, arrow_length_ratio=0.5)
# ax[:set_xlim](-2,2)
# ax[:set_ylim](-2,2)
# ax[:set_zlim](-2,2)
# ax[:set_xlabel]("x axis")
# ax[:set_ylabel]("y axis")
# ax[:set_zlabel]("z axis")
# fig[:show]()
