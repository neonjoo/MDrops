using CSV
using FastGaussQuadrature
using LinearAlgebra
using StatsBase
cd("/home/andris/MDrops/")
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")


points_csv= CSV.read("./meshes/points_ellipse_fewN.csv", header=0)
faces_csv = CSV.read("./meshes/faces_ellipse_fewN.csv", header=0)
println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')


edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)

(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)
k1,k2 = make_pc(CDE)
divn = k1+k2


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

        return zeta*r1 + xi*r2 + eta*r3
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

function make_izlgaussianF(points, faces, normals; gaussorder = 3)

    F = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number

            function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                n = B * zeta_xi_eta
                return n/norm(n)
            end

            function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = r/norm(r)^3 * (1 - 3*dot(r,nx)^2 / norm(r)^2)
                return f/(8pi)
            end

            function make_f_times_normr(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = r/norm(r)^2 * (1 - 3*dot(r,nx)^2 / norm(r)^2)
                return f/(8pi)
            end


            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                # singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])
                #
                # x1 = points[:,faces[singul_ind,i]]
                # x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # n1 = normals[:,faces[singul_ind,i]]
                # n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_izltrapezF(points, faces, normals)
    # replaces singular values with 0, should be same as F_old
    F = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number

            # function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
            #     nx = make_n(x,x1,x2,x3,n1,n2,n3)
            #     r = x - y
            #
            #     f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
            #     return f/(8pi)
            # end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]
                ny = normals[:,ykey]

                r1 = x1 - points[:,ykey]
                r2 = x2 - points[:,ykey]
                r3 = x3 - points[:,ykey]

                f1 = 1/(8pi) * ( dot(r1,n1)*ny + dot(r1,ny)*n1 + (1-dot(n1,ny))*r1 - 3*r1*dot(n1+ny,r1)*dot(r1,n1)/norm(r1)^2 )/ norm(r1)^3
                f2 = 1/(8pi) * ( dot(r2,n2)*ny + dot(r2,ny)*n2 + (1-dot(n2,ny))*r2 - 3*r2*dot(n2+ny,r2)*dot(r2,n2)/norm(r2)^2 )/ norm(r2)^3
                f3 = 1/(8pi) * ( dot(r3,n3)*ny + dot(r3,ny)*n3 + (1-dot(n3,ny))*r3 - 3*r3*dot(n3+ny,r3)*dot(r3,n3)/norm(r3)^2 )/ norm(r3)^3

                dS = 0.5 * norm(cross(x2-x1,x3-x1))

                F[:,ykey] += (f1+f2+f3)/3 * dS
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                # singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])
                #
                # x1 = points[:,faces[singul_ind,i]]
                # x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # n1 = normals[:,faces[singul_ind,i]]
                # n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # ny = normals[:,ykey]
                #
                # r1 = x1 - points[:,ykey]
                # r2 = x2 - points[:,ykey]
                # r3 = x3 - points[:,ykey]
                #
                # f1 = [0.,0.,0.]
                # f2 = 1/(8pi) * ( dot(r2,n2)*ny + dot(r2,ny)*n2 + (1-dot(n2,ny))*r2 - 3*r2*dot(n2+ny,r2)*dot(r2,n2)/norm(r2)^2 )/ norm(r2)^3
                # f3 = 1/(8pi) * ( dot(r3,n3)*ny + dot(r3,ny)*n3 + (1-dot(n3,ny))*r3 - 3*r3*dot(n3+ny,r3)*dot(r3,n3)/norm(r3)^2 )/ norm(r3)^3
                #
                # dS = 0.5 * norm(cross(x2-x1,x3-x1))
                #
                # F[:,ykey] += (f1+f2+f3)/3 * dS
            end

        end
    end
    return F
end

# g = make_izlgaussianF(points, faces, normals; gaussorder = 3)
# g10 = make_izlgaussianF(points, faces, normals; gaussorder = 10)
# t = make_izltrapezF(points, faces, normals)

function dbug(points, faces, normals; gaussorder = 3)

    F = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number

            function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                n = B * zeta_xi_eta
                return n/norm(n)
            end

            function make_lininterp(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters
                return B * zeta_xi_eta
            end

            function make_f(x,x1,x2,x3,n1,n2,n3,points,normals,ykey)
                y=points[:,ykey]
                ny=normals[:,ykey]
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                return 1/(8pi) * ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
            end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

                ny = normals[:,ykey]

                r1 = x1 - points[:,ykey]
                r2 = x2 - points[:,ykey]
                r3 = x3 - points[:,ykey]

                f1 = 1/(8pi) * ( dot(r1,n1)*ny + dot(r1,ny)*n1 + (1-dot(n1,ny))*r1 - 3*r1*dot(n1+ny,r1)*dot(r1,n1)/norm(r1)^2 )/ norm(r1)^3
                f2 = 1/(8pi) * ( dot(r2,n2)*ny + dot(r2,ny)*n2 + (1-dot(n2,ny))*r2 - 3*r2*dot(n2+ny,r2)*dot(r2,n2)/norm(r2)^2 )/ norm(r2)^3
                f3 = 1/(8pi) * ( dot(r3,n3)*ny + dot(r3,ny)*n3 + (1-dot(n3,ny))*r3 - 3*r3*dot(n3+ny,r3)*dot(r3,n3)/norm(r3)^2 )/ norm(r3)^3

                dS = 0.5 * norm(cross(x2-x1,x3-x1))

                gs = gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,points,normals,ykey), x1,x2,x3,gaussorder)
                lin_gs = gauss_nonsingular(x->make_lininterp(x,x1,x2,x3,f1,f2,f3), x1,x2,x3,gaussorder)
                tr = (f1+f2+f3)/3 * dS
                println("-------------------------")
                println("triangle number ",i)
                println("gaus   value ",gs)
                println("trap   value ",tr)
                println("lin_gs value ",lin_gs)
                readline()

            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                # singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])
                #
                # x1 = points[:,faces[singul_ind,i]]
                # x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # n1 = normals[:,faces[singul_ind,i]]
                # n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function dbugcrv(points, faces, normals,divn_p; gaussorder = 3)

    F = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number

            function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                n = B * zeta_xi_eta
                return n/norm(n)
            end

            function make_divn(x,x1,x2,x3,divn1,divn2,divn3) # divn linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [divn1 divn2 divn3] # vector of divn

                zeta_xi_eta = A \ x # find local triangle parameters

                return dot(B, zeta_xi_eta)
            end

            function make_lininterp(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters
                return B * zeta_xi_eta
            end

            function make_f(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3,points,ykey)
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                divn = make_divn(x,x1,x2,x3,divn1,divn2,divn3)
                y=points[:,ykey]
                r = x - y
                Gji = (Matrix(1.0I, (3,3)) / norm(r) +
                r * r' / norm(r)^3)

                return -( divn*Gji*nx)/(8pi)
            end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

                divn1 = divn_p[faces[1,i]]
                divn2 = divn_p[faces[2,i]]
                divn3 = divn_p[faces[3,i]]

                r1 = x1 - points[:,ykey]
                r2 = x2 - points[:,ykey]
                r3 = x3 - points[:,ykey]


                G1 = Matrix(1.0I, (3,3)) / norm(r1) + r1 * r1' / norm(r1)^3
                G2 = Matrix(1.0I, (3,3)) / norm(r2) + r2 * r2' / norm(r2)^3
                G3 = Matrix(1.0I, (3,3)) / norm(r3) + r3 * r3' / norm(r3)^3

                f1 = -( divn1*G1*n1)/(8pi)
                f2 = -( divn2*G2*n2)/(8pi)
                f3 = -( divn3*G3*n3)/(8pi)

                dS = 0.5 * norm(cross(x2-x1,x3-x1))

                gs = gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3,points,ykey), x1,x2,x3,gaussorder)
                lin_gs = gauss_nonsingular(x->make_lininterp(x,x1,x2,x3,f1,f2,f3), x1,x2,x3,gaussorder)
                tr = (f1+f2+f3)/3 * dS
                println("-------------------------")
                println("triangle number ",i)
                println("gaus   value ",gs)
                println("trap   value ",tr)
                println("lin_gs value ",lin_gs)
                readline()

            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                # singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])
                #
                # x1 = points[:,faces[singul_ind,i]]
                # x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # n1 = normals[:,faces[singul_ind,i]]
                # n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                # n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]
                #
                # F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

#F = dbugcrv(points, faces, normals,divn; gaussorder = 3)
F10 = dbugcrv(points, faces, normals,divn; gaussorder = 10)
