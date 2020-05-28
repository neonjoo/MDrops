using CSV
using FastGaussQuadrature
using LinearAlgebra
using StatsBase
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



points_csv= CSV.read("./meshes/points_ellipse_fewN.csv", header=0)
faces_csv = CSV.read("./meshes/faces_ellipse_fewN.csv", header=0)
println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')

points, faces = expand_icosamesh(R=1, depth=2)
points = Array{Float64}(points)
faces = Array{Int64}(faces)

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)

(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)
k1,k2 = make_pc(CDE)
divn = k1+k2

function make_gaussianF(points, faces, normals; gaussorder = 3)

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

            function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
                return f/(8pi)
            end

            function make_f_times_normr(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^2
                return f/(8pi)
            end


            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_nonsnggaussianF(points, faces, normals; gaussorder = 3)

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

            function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
                return f/(8pi)
            end

            function make_f_times_normr(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^2
                return f/(8pi)
            end


            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                #F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_snggaussianF(points, faces, normals; gaussorder = 3)

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

            function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = r/norm(r)^3 * (1 - 3*dot(r,nx)^2 / norm(r)^2)
                return f/(8pi)
            end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                F[:,ykey] += gauss_weaksingular(x->norm(x-x1)*make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
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

            function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
                return f/(8pi)
            end

            function make_f_times_normr(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^2
                return f/(8pi)
            end


            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                #F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_curvgaussianF(points, faces, normals, divn_p; gaussorder = 3)

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

                divn = dot(B, zeta_xi_eta)
                return divn
            end

            function make_f(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3; y=points[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                divn = make_divn(x,x1,x2,x3,divn1,divn2,divn3)
                r = x - y
                Gji = (Matrix(1.0I, (3,3)) / norm(r) +
                r * r' / norm(r)^3)

                f = -( divn*Gji*nx)
                return f/(8pi)
            end

            function make_f_times_normr(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3; y=points[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                divn = make_divn(x,x1,x2,x3,divn1,divn2,divn3)
                r = x - y
                Gji_times_normr = (Matrix(1.0I, (3,3))  +
                r * r' / norm(r)^2)

                f = -( divn*Gji_times_normr*nx)
                return f/(8pi)
            end


            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                divn1 = divn_p[faces[1,i]]
                divn2 = divn_p[faces[2,i]]
                divn3 = divn_p[faces[3,i]]


                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                divn1 = divn_p[faces[singul_ind,i]]
                divn2 = divn_p[faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                divn3 = divn_p[faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_crvlssgaussianF(points, faces, normals; gaussorder = 3)

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
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_specialgaussianF(points, faces, normals; gaussorder = 3)

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

            function make_f(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ( dot(r,nx)*ny + dot(r,ny)*nx + (1-dot(nx,ny))*r - 3*r*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
                return f/(8pi)
            end

            function make_fny(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                # component only in ny direction, becaus in the limit x->y,
                # the function is antisymmetric in the tangential direction
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ny*( dot(r,nx) + dot(r,ny)*dot(nx,ny) + (1-dot(nx,ny))*dot(r,ny) - 3*dot(r,ny)*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^3
                return f/(8pi)
            end

            function make_fny_timesnormr(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
                # component only in ny direction, becaus in the limit x->y,
                # the function is antisymmetric in the tangential direction
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = x - y

                f = ny*( dot(r,nx) + dot(r,ny)*dot(nx,ny) + (1-dot(nx,ny))*dot(r,ny) - 3*dot(r,ny)*dot(nx+ny,r)*dot(r,nx)/norm(r)^2 )/ norm(r)^2
                return f/(8pi)
            end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1 = normals[:,faces[1,i]]
                n2 = normals[:,faces[2,i]]
                n3 = normals[:,faces[3,i]]

                F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            else # if is singular triangle

                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                F[:,ykey] += gauss_weaksingular(x->make_fny_timesnormr(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_fny(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

specialF = make_specialgaussianF(points, faces, normals, gaussorder = 10)
specialFweak = make_specialgaussianF(points, faces, normals, gaussorder = 10)
