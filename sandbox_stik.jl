using CSV
using FastGaussQuadrature
using LinearAlgebra
using StatsBase
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")



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

#points, faces = expand_icosamesh(R=1, depth=2)
#points = Array{Float64}(points)
#faces = Array{Int64}(faces)

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

function solidbody_project(points, faces, w)
    # projects the surfaces velocity w
    # on the rotation and translations of a solid body

    function make_w(x,x1,x2,x3,w1,w2,w3) # vector linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [w1 w2 w3] # matrix of vertex the vectorfunction

        zeta_xi_eta = A \ x # find local triangle parameters

        return B * zeta_xi_eta
    end

    S = 0
    yc = [0.,0.,0.]
    V = [0.,0.,0.]
    Omega = [0.,0.,0.]
    M = zeros(3,3)

    for i = 1:size(faces,2) # triangle number
        x1 = points[:,faces[1,i]]
        x2 = points[:,faces[2,i]]
        x3 = points[:,faces[3,i]]

        w1 = w[:,faces[1,i]]
        w2 = w[:,faces[2,i]]
        w3 = w[:,faces[3,i]]

        deltaS = norm(cross(x2-x1,x3-x1))/2
        S += deltaS
        # trapezoid rule because linear functions
        yc += (x1+x2+x3)/3 * deltaS
        V += (w1+w2+w3)/3 * deltaS
    end

    yc /= S
    V /= S

    function x_tilde_cross_w(x,yc,w1,w2,w3,x1,x2,x3)
        intw = make_w(x,x1,x2,x3,w1,w2,w3)
        intxtilde = x - yc

        return cross(intxtilde,intw)
    end

    function Mfun(x,yc)
        intxtilde = x - yc
        return dot(intxtilde,intxtilde) * Matrix(1.0I,3,3) - intxtilde * intxtilde'
    end

    for i = 1:size(faces,2) # triangle number
        x1 = points[:,faces[1,i]]
        x2 = points[:,faces[2,i]]
        x3 = points[:,faces[3,i]]

        w1 = w[:,faces[1,i]]
        w2 = w[:,faces[2,i]]
        w3 = w[:,faces[3,i]]

        deltaS = norm(cross(x2-x1,x3-x1))/2

        # gaussian integration, because quadratic

        Omega += gauss_nonsingular(x->x_tilde_cross_w(x,yc,w1,w2,w3,x1,x2,x3),x1,x2,x3,2)
        M += gauss_nonsingular(x->Mfun(x,yc),x1,x2,x3,2)
    end

    Omega = inv(M) * Omega

    wprim = Array{Float64}(undef, size(points))
    for ykey = 1:size(points,2)
        wprim[:,ykey] = V + cross(Omega,points[:,ykey]-yc)
    end

    return wprim
end

function make_wielandtL(points, faces, normals, w, lambda; gaussorder = 3)
    # w = Lw + F


    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_w(x,x1,x2,x3,w1,w2,w3) # vector linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [w1 w2 w3] # matrix of vertex the vectorfunction

        zeta_xi_eta = A \ x # find local triangle parameters

        return B * zeta_xi_eta
    end

    function make_f(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey)
        # Tijk = -6 rrr/|r|^5
        y=points[:,ykey]
        ny=normals[:,ykey]
        wy = w[:,ykey]

        nx = make_n(x,x1,x2,x3,n1,n2,n3)
        wx = make_w(x,x1,x2,x3,w1,w2,w3) # just the same interpolation for w

        r = x - y
        return 1/(4pi)*(-6)*dot(wx-wy,r)*r*dot(r,nx)/norm(r)^5 - ny/S*dot(wx,nx)
    end

    S = make_S(points,faces)

    L = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number
            #println(i)
            #println(ykey)
            x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
            #println("xok")

            #println(size(w))
            w1, w2, w3 = [w[:, faces[j,i]] for j in 1:3]
            #println("wok")

            n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
            #println("nok")
            #println("fails above")

            hold = gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey), x1,x2,x3,gaussorder)
            #println("deet")
            #println(hold)
            #println(L[:,ykey])
            L[:,ykey] += hold
            #println("endisnigh")
        end
    end

    L += solidbody_project(points, faces, w)

    return (1-lambda)/2 * L
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

function fastsolidbody_project(points, faces, w)
    # projects the surfaces velocity w
    # on the rotation and translations of a solid body
    # w = 0 0 0 0 0 1 0 0 , the position of 1 is encoded by oneind
    oneind = findfirst(x->x==1,w)[1]

    function make_w(x,x1,x2,x3,w1,w2,w3) # vector linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [w1 w2 w3] # matrix of vertex the vectorfunction

        zeta_xi_eta = A \ x # find local triangle parameters

        return B * zeta_xi_eta
    end

    S = 0
    yc = [0.,0.,0.]
    V = [0.,0.,0.]
    Omega = [0.,0.,0.]
    M = zeros(3,3)

    for i = 1:size(faces,2) # triangle number

        x1 = points[:,faces[1,i]]
        x2 = points[:,faces[2,i]]
        x3 = points[:,faces[3,i]]

        deltaS = norm(cross(x2-x1,x3-x1))/2

        if oneind in faces[:,i]
            w1 = w[:,faces[1,i]]
            w2 = w[:,faces[2,i]]
            w3 = w[:,faces[3,i]]
            # trapezoid rule because linear functions
            V += (w1+w2+w3)/3 * deltaS
        end


        S += deltaS
        # trapezoid rule because linear functions
        yc += (x1+x2+x3)/3 * deltaS
    end

    yc /= S
    V /= S

    function x_tilde_cross_w(x,yc,w1,w2,w3,x1,x2,x3)
        intw = make_w(x,x1,x2,x3,w1,w2,w3)
        intxtilde = x - yc

        return cross(intxtilde,intw)
    end

    function Mfun(x,yc)
        intxtilde = x - yc
        return dot(intxtilde,intxtilde) * Matrix(1.0I,3,3) - intxtilde * intxtilde'
    end

    for i = 1:size(faces,2) # triangle number
        x1 = points[:,faces[1,i]]
        x2 = points[:,faces[2,i]]
        x3 = points[:,faces[3,i]]

        deltaS = norm(cross(x2-x1,x3-x1))/2

        if oneind in faces[:,i]
            w1 = w[:,faces[1,i]]
            w2 = w[:,faces[2,i]]
            w3 = w[:,faces[3,i]]
            # gaussian integration, because quadratic
            Omega += gauss_nonsingular(x->x_tilde_cross_w(x,yc,w1,w2,w3,x1,x2,x3),x1,x2,x3,2)
        end

        # gaussian integration, because quadratic
        M += gauss_nonsingular(x->Mfun(x,yc),x1,x2,x3,2)
    end

    Omega = inv(M) * Omega

    wprim = Array{Float64}(undef, size(points))
    for ykey = 1:size(points,2)
        wprim[:,ykey] = V + cross(Omega,points[:,ykey]-yc)
    end

    return wprim
end

function make_fastwielandtL(points, faces, normals, w, lambda; gaussorder = 3)
    # w = Lw + F
    # w = 0 0 0 0 0 1 0 0 , the position of 1 is encoded by oneind
    oneind = findfirst(x->x==1,w)[1]

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_w(x,x1,x2,x3,w1,w2,w3) # vector linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [w1 w2 w3] # matrix of vertex the vectorfunction

        zeta_xi_eta = A \ x # find local triangle parameters

        return B * zeta_xi_eta
    end

    function make_f(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey)
        # Tijk = -6 rrr/|r|^5
        y=points[:,ykey]
        ny=normals[:,ykey]
        wy = w[:,ykey]

        nx = make_n(x,x1,x2,x3,n1,n2,n3)
        wx = make_w(x,x1,x2,x3,w1,w2,w3) # just the same interpolation for w

        r = x - y
        return 1/(4pi)*(-6)*dot(wx-wy,r)*r*dot(r,nx)/norm(r)^5 - ny/S*dot(wx,nx)
    end

    S = make_S(points,faces)

    L = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number
            if oneind in faces[:,i]
                #println(i)
                #println(ykey)
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                #println("xok")

                #println(size(w))
                w1, w2, w3 = [w[:, faces[j,i]] for j in 1:3]
                #println("wok")

                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
                #println("nok")
                #println("fails above")

                hold = gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey), x1,x2,x3,gaussorder)
                #println("deet")
                #println(hold)
                #println(L[:,ykey])
                L[:,ykey] += hold
                #println("endisnigh")
            end
        end
    end

    L += fastsolidbody_project(points, faces, w)

    return (1-lambda)/2 * L
end

function make_fastwielandtLsng(points, faces, normals, w, lambda; gaussorder = 3)
    # w = Lw + F
    # w = 0 0 0 0 0 1 0 0 , the position of 1 is encoded by oneind
    oneind = findfirst(x->x==1,w)[1]

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_w(x,x1,x2,x3,w1,w2,w3) # vector linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [w1 w2 w3] # matrix of vertex the vectorfunction

        zeta_xi_eta = A \ x # find local triangle parameters

        return B * zeta_xi_eta
    end

    function make_f(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey)
        # Tijk = -6 rrr/|r|^5
        y=points[:,ykey]
        ny=normals[:,ykey]
        wy = w[:,ykey]

        nx = make_n(x,x1,x2,x3,n1,n2,n3)
        wx = make_w(x,x1,x2,x3,w1,w2,w3) # just the same interpolation for w

        r = x - y
        return 1/(4pi)*(-6)*dot(wx-wy,r)*r*dot(r,nx)/norm(r)^5 - ny/S*dot(wx,nx)
    end

    function make_ftimesnormr(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey)
        # Tijk = -6 rrr/|r|^5
        y=points[:,ykey]
        ny=normals[:,ykey]
        wy = w[:,ykey]

        nx = make_n(x,x1,x2,x3,n1,n2,n3)
        wx = make_w(x,x1,x2,x3,w1,w2,w3) # just the same interpolation for w

        r = x - y
        return 1/(4pi)*(-6)*dot(wx-wy,r)*r*dot(r,nx)/norm(r)^4 - ny/S*dot(wx,nx)*norm(r)
    end

    S = make_S(points,faces)

    L = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number
            if oneind in faces[:,i]
                if !(ykey in faces[:,i]) # if not singular triangle


                    x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                    w1, w2, w3 = [w[:, faces[j,i]] for j in 1:3]
                    n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

                    L[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey), x1,x2,x3,gaussorder)

                else
                    singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                    x1 = points[:,faces[singul_ind,i]]
                    x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                    x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                    w1 = w[:,faces[singul_ind,i]]
                    w2 = w[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                    w3 = w[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                    n1 = normals[:,faces[singul_ind,i]]
                    n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                    n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                    L[:,ykey] += gauss_weaksingular(x->make_ftimesnormr(x,x1,x2,x3,n1,n2,n3,w1,w2,w3,S,ykey), x1,x2,x3,gaussorder)
                end
            end
        end
    end

    L += fastsolidbody_project(points, faces, w)

    return (1-lambda)/2 * L
end

function make_zerotrapezF(points, faces, normals)
    # replaces singular values with 0, should be same as F_old
    F = zeros(size(points))
    for ykey = 1:size(points,2)
        for i = 1:size(faces,2) # triangle number

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
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind + 1 - 1) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 2 - 1) % 3 + 1,i]]

                ny = normals[:,ykey]

                r1 = x1 - points[:,ykey]
                r2 = x2 - points[:,ykey]
                r3 = x3 - points[:,ykey]

                f1 = [0.,0.,0.]
                f2 = 1/(8pi) * ( dot(r2,n2)*ny + dot(r2,ny)*n2 + (1-dot(n2,ny))*r2 - 3*r2*dot(n2+ny,r2)*dot(r2,n2)/norm(r2)^2 )/ norm(r2)^3
                f3 = 1/(8pi) * ( dot(r3,n3)*ny + dot(r3,ny)*n3 + (1-dot(n3,ny))*r3 - 3*r3*dot(n3+ny,r3)*dot(r3,n3)/norm(r3)^2 )/ norm(r3)^3

                dS = 0.5 * norm(cross(x2-x1,x3-x1))

                F[:,ykey] += (f1+f2+f3)/3 * dS
            end

        end
    end
    return F
end

function make_izltrapezF(points, faces, normals)
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

function make_zerogaussianF(points, faces, normals; gaussorder = 3)

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

                ny = normals[:,ykey]

                r1 = x1 - points[:,ykey]
                r2 = x2 - points[:,ykey]
                r3 = x3 - points[:,ykey]

                f1 = [0.,0.,0.]
                f2 = 1/(8pi) * ( dot(r2,n2)*ny + dot(r2,ny)*n2 + (1-dot(n2,ny))*r2 - 3*r2*dot(n2+ny,r2)*dot(r2,n2)/norm(r2)^2 )/ norm(r2)^3
                f3 = 1/(8pi) * ( dot(r3,n3)*ny + dot(r3,ny)*n3 + (1-dot(n3,ny))*r3 - 3*r3*dot(n3+ny,r3)*dot(r3,n3)/norm(r3)^2 )/ norm(r3)^3

                dS = 0.5 * norm(cross(x2-x1,x3-x1))

                F[:,ykey] += (f1+f2+f3)/3 * dS
            end

        end
    end
    return F
end

function make_justsingF(points, faces, normals; gaussorder = 3)

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

            if !(ykey in faces[:,i]) # if not singular triangle
                # x1 = points[:,faces[1,i]]
                # x2 = points[:,faces[2,i]]
                # x3 = points[:,faces[3,i]]
                #
                # n1 = normals[:,faces[1,i]]
                # n2 = normals[:,faces[2,i]]
                # n3 = normals[:,faces[3,i]]
                #
                # F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
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

                ny = normals[:,ykey]

                r1 = x1 - points[:,ykey]
                r2 = x2 - points[:,ykey]
                r3 = x3 - points[:,ykey]

                f1 = [0.,0.,0.]
                f2 = 1/(8pi) * ( dot(r2,n2)*ny + dot(r2,ny)*n2 + (1-dot(n2,ny))*r2 - 3*r2*dot(n2+ny,r2)*dot(r2,n2)/norm(r2)^2 )/ norm(r2)^3
                f3 = 1/(8pi) * ( dot(r3,n3)*ny + dot(r3,ny)*n3 + (1-dot(n3,ny))*r3 - 3*r3*dot(n3+ny,r3)*dot(r3,n3)/norm(r3)^2 )/ norm(r3)^3

                dS = 0.5 * norm(cross(x2-x1,x3-x1))

                F[:,ykey] += (f1+f2+f3)/3 * dS
            end

        end
    end
    return F
end

function make_justgaussingF(points, faces, normals, divn_p; gaussorder = 3)

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


                # F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3), x1,x2,x3,gaussorder)
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

function make_justgausnonsingF(points, faces, normals, divn_p; gaussorder = 3)

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

                #F[:,ykey] += gauss_weaksingular(x->make_f_times_normr(x,x1,x2,x3,n1,n2,n3,divn1,divn2,divn3), x1,x2,x3,gaussorder)
                #F[:,ykey] += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
            end

        end
    end
    return F
end

function make_justnonsingF(points, faces, normals; gaussorder = 3)

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
# g = 7
# order = 1
# println(gauss_nonsingular(x->dot(x,x)^(order/2), [0.,0.,0.],[0.,0.,1.],[0., 5., 0.5],g))
# print((gauss_nonsingular(x->dot(x,x)^(order/2), [0.,0.,0.],[0.,0.,1.],[0., 5., 0.5],g)
# - gauss_nonsingular(x->dot(x,x)^(order/2), [0.,0.,0.],[0.,0.,1.],[0., 5., 0.5],10))
# / gauss_nonsingular(x->dot(x,x)^(order/2), [0.,0.,0.],[0.,0.,1.],[0., 5., 0.5],10)*100)
# println(" %")
# println()
# r1 = [0.,0.,0.]
# r2 = [0.,0.,1.]
# r3 = [0., 1., 0.]
# println((trapezoid_nonsingular(dot(r1,r1),dot(r2,r2),dot(r3,r3),r1,r2,r3) - gauss_nonsingular(x->dot(x,x), r1,r2,r3,10)) / gauss_nonsingular(x->dot(x,x), r1,r2,r3,10) *100 )
#
#
# function trapezoid_nonsingular(f1,f2,f3,r1,r2,r3)
#     #f - function of r
#     deltaS = norm(cross(r2-r1,r3-r1))/2
#     return (f1+f2+f3)/3 * deltaS
# end
#
# println(trapezoid_nonsingular(r1,r2,r3,r1,r2,r3))
# println(gauss_nonsingular(x->x, r1,r2,r3,2))
# println()


# w = Lw + F
# L = LI
# 0 = (L-I)w + F
# w = -(L-I) \ F
# v = w + (lambda-1)/2 * solidbody_project(w)
#velocities = reshape(varr ,3,size(vertices,2))
#varr = reshape(velocities ,1,3*size(vertices,2))
lambda = 10
Lmat  = Array{Float64}(undef, (3*size(points,2),3*size(points,2)))
for i = 1:3*size(points,2)
    varr = zeros(1,3*size(points,2))
    varr[i] = 1.
    vels = reshape(varr ,3,size(points,2))
    vels = make_fastwielandtLsng(points, faces, normals, vels, lambda; gaussorder = 3)
    varr = reshape(vels ,1,3*size(points,2))
    Lmat[:,i] = varr
end

F = make_curvgaussianF(points, faces, normals, divn; gaussorder = 3)
F_izl = make_izlgaussianF(points, faces, normals; gaussorder = 3)
F10 = make_curvgaussianF(points, faces, normals, divn; gaussorder = 10)
Farr = reshape(F ,1,3*size(points,2))
warr = -(Lmat-(Matrix(1.0I, size(Lmat)) )) \ Farr'
w = reshape(warr ,3,size(points,2))
velocities_sng = w + (lambda-1)/2 * solidbody_project(points,faces,w)

velocities_test = make_magvelocities(points, normals, lambda, 0, 0, zeros(size(points)), zeros(size(points)))
F_test = make_magvelocities(points, normals, 1, 0, 0, zeros(size(points)), zeros(size(points)))


# check if divergence is zero
deltaS = make_dS(points,faces)
V = 1/3 * sum(sum(points .* normals,dims=1) .* deltaS')
Fcheck = sum(sum(F .* normals,dims=1) .* deltaS') / V
F10check = sum(sum(F10 .* normals,dims=1) .* deltaS') / V
F_testcheck = sum(sum(F_test .* normals,dims=1) .* deltaS') / V
velocities_sngcheck = sum(sum(velocities_sng .* normals,dims=1) .* deltaS') / V
velocities_testcheck = sum(sum(velocities_test .* normals,dims=1) .* deltaS') / V


#
# using PyPlot
# pygui()
#
# fig = figure(figsize=(7,7))
# ax = fig[:gca](projection="3d")
#
# (x, y, z) = [points[i,:] for i in 1:3]
# (vx, vy, vz) = [velocities_sng[i,:] for i in 1:3]
#
# ax[:scatter](x,y,z, s=2,color="k")
# ax[:quiver](x,y,z,vx,vy,vz, length=30, arrow_length_ratio=0.5)
#
# ax[:set_xlim](-2,2)
# ax[:set_ylim](-2,2)
# ax[:set_zlim](-2,2)
# ax[:set_xlabel]("x axis")
# ax[:set_ylabel]("y axis")
# ax[:set_zlabel]("z axis")
# fig[:show]()
#
