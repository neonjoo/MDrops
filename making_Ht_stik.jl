using CSV
using FastGaussQuadrature
using LinearAlgebra
using StatsBase
using QuadGK
using Cubature
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

function hquad_weaksingular_vec(q::Function,r1,r2,r3)
    # integral of q(r)/|r-r1|, where q(r) is nonsingular
    # singularity on r1

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

    function make_integrand(uv, B, C, r1,r2,r3)
        chi = make_chi(uv[1])
        make_R(uv[1]) / sqrt( cos(chi)^2 + B*sin(2*chi) + C*sin(chi)^2 ) * q( make_r(uv[1],uv[2],r1,r2,r3) )
    end

    intdim = 3 #output is a vector#
    #println(make_integrand([0.5, 0.5], B, C, r1,r2,r3))
    v = [0, 0, 0]
    (intval, err) = hcubature(intdim, (uv,v) -> v[:] = make_integrand(uv, B, C, r1,r2,r3), [-1, -1], [1, 1];
                      reltol=1e-8, abstol=1e-8, maxevals=0,
                      error_norm = Cubature.INDIVIDUAL)
    println("intval")
    println(intval)
    println("error in h quadrature:")
    println(err)
    println()
    #readline()
    intval *= pi*hs / (8 * norm(r2-r1))
    return intval
end

function gauss_flat_local_vec(f::Function,r1,r2,r3,gaussorder)
    # integral of f(r)
    # singularity on r1

    n_triang = cross(r3-r2,r3-r1) / norm(cross(r3-r2,r3-r1))

    r1_loc, r2_loc, r3_loc = r1-r1, r2-r1, r3-r1
    r2_loc, r3_loc = to_local(r2_loc,n_triang), to_local(r3_loc,n_triang)

    x2, x3 = r2_loc[1], r3_loc[1]
    y2, y3 = r2_loc[2], r3_loc[2]

    # angles are from 0 to 2pi, masured from the x axis
    theta2, theta3 = mod2pi(atan(y2,x2)+2pi + 1e-15) - 1e-15, mod2pi(atan(y3,x3)+2pi + 1e-15) - 1e-15

    u, wu = gausslegendre(gaussorder) # u from -1 to 1
    v, wv = u, wu # for convenient notation

    function make_theta(v, theta2, theta3)
        return (v + 1)/2 * (theta3 - theta2) + theta2
    end

    function make_rho_m(theta, x2, x3, y2, y3)
        return (y2 - x2 * (y3-y2)/(x3-x2)) / (sin(theta) - cos(theta) * (y3-y2)/(x3-x2))
    end

    function make_rho(u, rho_m)
        return (u + 1)/2 * rho_m
    end

    function make_r(rho,theta, n_triang, r1)
        x = rho * cos(theta)
        y = rho * sin(theta)
        z = 0.

        r_loc = [x, y, z]

        return to_global(r_loc,n_triang) + r1
    end

    intval = [0.,0.,0.]

    for i = 1:length(v)
        theta = make_theta(v[i], theta2, theta3)
        rho_m = make_rho_m(theta, x2, x3, y2, y3)
        for k = 1:length(u)
            rho = make_rho(u[k], rho_m)
            intval += wv[i]*(theta3 - theta2)/2 *
                    wu[k]*f( make_r(rho,theta, n_triang, r1) ) * rho * rho_m/2

        end
    end

    if theta2>theta3 # change of sign if integration limits are reversed
        intval *= -1
    end

    return intval
end

function gauss_flat_local_scal(f::Function,r1,r2,r3,gaussorder)
    # integral of f(r)
    # singularity on r1

    n_triang = cross(r3-r2,r3-r1) / norm(cross(r3-r2,r3-r1))

    r1_loc, r2_loc, r3_loc = r1-r1, r2-r1, r3-r1
    r2_loc, r3_loc = to_local(r2_loc,n_triang), to_local(r3_loc,n_triang)

    x2, x3 = r2_loc[1], r3_loc[1]
    y2, y3 = r2_loc[2], r3_loc[2]

    # angles are from 0 to 2pi, masured from the x axis
    theta2, theta3 = mod2pi(atan(y2,x2)+2pi + 1e-15) - 1e-15, mod2pi(atan(y3,x3)+2pi + 1e-15) - 1e-15

    u, wu = gausslegendre(gaussorder) # u from -1 to 1
    v, wv = u, wu # for convenient notation

    function make_theta(v, theta2, theta3)
        return (v + 1)/2 * (theta3 - theta2) + theta2
    end

    function make_rho_m(theta, x2, x3, y2, y3)
        return (y2 - x2 * (y3-y2)/(x3-x2)) / (sin(theta) - cos(theta) * (y3-y2)/(x3-x2))
    end

    function make_rho(u, rho_m)
        return (u + 1)/2 * rho_m
    end

    function make_r(rho,theta, n_triang, r1)
        x = rho * cos(theta)
        y = rho * sin(theta)
        z = 0.

        r_loc = [x, y, z]

        return to_global(r_loc,n_triang) + r1
    end

    intval = 0.

    for i = 1:length(v)
        theta = make_theta(v[i], theta2, theta3)
        rho_m = make_rho_m(theta, x2, x3, y2, y3)
        for k = 1:length(u)
            rho = make_rho(u[k], rho_m)
            intval += wv[i]*(theta3 - theta2)/2 *
                    wu[k]*f( make_r(rho,theta, n_triang, r1) ) * rho * rho_m/2

        end
    end

    if theta2>theta3 # change of sign if integration limits are reversed
        intval *= -1
    end

    return intval
end

function make_L_sing(points, faces, normals; gaussorder=5)
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

function make_deltaH_normal(points, faces, normals, mu, H0)
    alpha = mu
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = make_L_sing(points, faces, normals; gaussorder=5)
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

    return deltaH_normal'
end

function make_H_tang(points, faces, normals, delta_H_normal, H0; gaussorder = 5)
    N = size(points, 2)
    H_tang = Array{Float64}(undef, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta_H_normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [dHn1 dHn2 dHn3] # vector of vertex delta_H_normal

        zeta_xi_eta = A \ x # find local triangle parameters

        return dot(B, zeta_xi_eta)
    end

    for ykey in 1:N
        n_cross_H = cross(normals[:, ykey], H0)

        for i in 1:size(faces, 2)

            function make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^3 * cross(dHnx*ny - dHny*nx , r)
            end

            function make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey], ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^2 * cross(dHnx*ny - dHny*nx , r)
            end

            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
                dHn1,dHn2,dHn3 = [delta_H_normal[faces[j,i]] for j in 1:3]

                n_cross_H += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 1) % 3 + 1,i]]

                dHn1 = delta_H_normal[faces[singul_ind,i]]
                dHn2 = delta_H_normal[faces[(singul_ind) % 3 + 1,i]]
                dHn3 = delta_H_normal[faces[(singul_ind + 1) % 3 + 1,i]]

                n_cross_H += gauss_weaksingular(x->make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
            end
        end
        H_tang[ykey] = norm(n_cross_H)
    end

    return H_tang'
end

function make_H_tang_doublecross(points, faces, normals, delta_H_normal, H0; gaussorder = 5)
    N = size(points, 2)
    H_tang = Array{Float64}(undef, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta_H_normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [dHn1 dHn2 dHn3] # vector of vertex delta_H_normal

        zeta_xi_eta = A \ x # find local triangle parameters

        return dot(B, zeta_xi_eta)
    end

    for ykey in 1:N
        n_cross_H = cross(cross(normals[:, ykey], H0), normals[:, ykey])

        for i in 1:size(faces, 2)

            function make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^3 * cross(cross(dHnx*ny - dHny*nx , r),ny)
            end

            function make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey], ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^2 * cross(cross(dHnx*ny - dHny*nx , r),ny)
            end

            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
                dHn1,dHn2,dHn3 = [delta_H_normal[faces[j,i]] for j in 1:3]

                n_cross_H += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 1) % 3 + 1,i]]

                dHn1 = delta_H_normal[faces[singul_ind,i]]
                dHn2 = delta_H_normal[faces[(singul_ind) % 3 + 1,i]]
                dHn3 = delta_H_normal[faces[(singul_ind + 1) % 3 + 1,i]]

                n_cross_H += gauss_weaksingular(x->make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
            end
        end
        H_tang[ykey] = norm(n_cross_H)
    end

    return H_tang'
end

function make_H_tang_semianal(points, faces, normals, delta_H_normal, H0, CDE; gaussorder = 5)
    N = size(points, 2)
    H_tang = Array{Float64}(undef, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta_H_normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [dHn1 dHn2 dHn3] # vector of vertex delta_H_normal

        zeta_xi_eta = A \ x # find local triangle parameters

        return dot(B, zeta_xi_eta)
    end

    for ykey in 1:N
        n_cross_H = cross(normals[:, ykey], H0)

        Rsq = zeros(3,3)
        C = zeros(3)
        rmean = 0.
        rmeancounter = 0
        for i in 1:size(faces, 2)

            function make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^3 * cross(dHnx*ny - dHny*nx , r)
            end
            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
                dHn1,dHn2,dHn3 = [delta_H_normal[faces[j,i]] for j in 1:3]

                n_cross_H += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                dHn1 = delta_H_normal[faces[singul_ind,i]]
                dHn2 = delta_H_normal[faces[(singul_ind) % 3 + 1,i]]
                dHn3 = delta_H_normal[faces[(singul_ind + 1) % 3 + 1,i]]

                r2 = to_local(x2 - x1, normals[:,ykey])
                r3 = to_local(x3 - x1, normals[:,ykey])

                rmean += norm(r2)
                rmean += norm(r3)

                rmeancounter += 2

                Rsq += 0.5 * r2 * r2'
                Rsq += 0.5 * r3 * r3'

                C += 0.5 * ( dHn2 * r2 - dHn1 * r2)
                C += 0.5 * ( dHn3 * r3 - dHn1 * r3)

            end
        end

        rmean /= rmeancounter
        Htilde = C \ Rsq
        #n_cross_H += pi/4 * to_global([Htilde[2] * rmean^4, -Htilde[1] * rmean^4, 0], normals[:,ykey])

        a, b, c = 2*CDE[1,ykey], CDE[2,ykey], 2*CDE[3,ykey]

        n_cross_H +=to_global(
        [Htilde[2]*rmean/4 - 1/256*(4b*(a+c)*Htilde[1] + (a^2 + 4b^2 + 2a*c + 5c^2)*Htilde[2])*rmean^3,
        -Htilde[1]*rmean/4 + 1/256*( (5a^2 + 4b^2 + 2a*c + c^2)*Htilde[1] + 4b*(a+c)*Htilde[2])*rmean^3,
         0]
        , normals[:,ykey])

        H_tang[ykey] = norm(n_cross_H)
    end

    return H_tang'
end

function make_H_tang_hquad(points, faces, normals, delta_H_normal, H0; gaussorder = 5)
    N = size(points, 2)
    H_tang = Array{Float64}(undef, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta_H_normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [dHn1 dHn2 dHn3] # vector of vertex delta_H_normal

        zeta_xi_eta = A \ x # find local triangle parameters

        return dot(B, zeta_xi_eta)
    end

    for ykey in 1:N
        n_cross_H = cross(normals[:, ykey], H0)

        for i in 1:size(faces, 2)

            function make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^3 * cross(dHnx*ny - dHny*nx , r)
            end

            function make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey], ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^2 * cross(dHnx*ny - dHny*nx , r)
            end

            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
                dHn1,dHn2,dHn3 = [delta_H_normal[faces[j,i]] for j in 1:3]

                n_cross_H += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 1) % 3 + 1,i]]

                dHn1 = delta_H_normal[faces[singul_ind,i]]
                dHn2 = delta_H_normal[faces[(singul_ind) % 3 + 1,i]]
                dHn3 = delta_H_normal[faces[(singul_ind + 1) % 3 + 1,i]]

                n_cross_H += hquad_weaksingular_vec(x->make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3)
            end
        end
        H_tang[ykey] = norm(n_cross_H)
    end

    return H_tang'
end

function make_H_tang_skip_singul(points, faces, normals, delta_H_normal, H0; gaussorder = 5)
    N = size(points, 2)
    H_tang = Array{Float64}(undef, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta_H_normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [dHn1 dHn2 dHn3] # vector of vertex delta_H_normal

        zeta_xi_eta = A \ x # find local triangle parameters

        return dot(B, zeta_xi_eta)
    end

    for ykey in 1:N
        n_cross_H = cross(normals[:, ykey], H0)

        for i in 1:size(faces, 2)

            function make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^3 * cross(dHnx*ny - dHny*nx , r)
            end

            function make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey], ny=normals[:,ykey], dHny = delta_H_normal[ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                r = y - x

                return 1/(4pi) * 1/norm(r)^2 * cross(dHnx*ny - dHny*nx , r)
            end

            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]
                dHn1,dHn2,dHn3 = [delta_H_normal[faces[j,i]] for j in 1:3]

                n_cross_H += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 1) % 3 + 1,i]]

                dHn1 = delta_H_normal[faces[singul_ind,i]]
                dHn2 = delta_H_normal[faces[(singul_ind) % 3 + 1,i]]
                dHn3 = delta_H_normal[faces[(singul_ind + 1) % 3 + 1,i]]

                #n_cross_H += hquad_weaksingular_vec(x->make_f_times_norm_r(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3)
            end
        end
        H_tang[ykey] = norm(n_cross_H)
    end

    return H_tang'
end

function demag_coefs(a, b, c)
    upper_limit = 2000
    Rq2(q) = (a^2+q) * (b^2+q) * (c^2+q)

    Nx = a*b*c/2 * quadgk(s -> 1/(a^2+s) / sqrt(Rq2(s)), 0, upper_limit)[1]
    Ny = a*b*c/2 * quadgk(s -> 1/(b^2+s) / sqrt(Rq2(s)), 0, upper_limit)[1]
    Nz = a*b*c/2 * quadgk(s -> 1/(c^2+s) / sqrt(Rq2(s)), 0, upper_limit)[1]

    return [Nx, Ny, Nz]
end

function field_theor(a, b, c, mu, H0)

    Ns = demag_coefs(a,b,c)

    Hx = H0[1] / (1 + Ns[1] * (mu-1))
    Hy = H0[2] / (1 + Ns[2] * (mu-1))
    Hz = H0[3] / (1 + Ns[3] * (mu-1))

    return [Hx, Hy, Hz]
end


points, faces = expand_icosamesh(R=1, depth=2)
points = Array{Float64}(points)
faces = Array{Int64}(faces)
normals = Normals(points, faces)
edges = make_edges(faces)
connectivity = make_connectivity(edges)
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)
H0 = [0., 0., 1.]
mu = 10
## test if normals are outward
if minimum(sum(normals .* points, dims=1)) > 0
    println("normals OK")
end


deltaHn = make_deltaH_normal(points, faces, normals, mu, H0)
Hn = deltaHn / (mu-1)
Ht = make_H_tang(points, faces, normals, deltaHn, H0; gaussorder = 5)
Ht_dc = make_H_tang_doublecross(points, faces, normals, deltaHn, H0; gaussorder = 5)
Ht_sa = make_H_tang_semianal(points, faces, normals, deltaHn, H0, CDE; gaussorder = 5)
Ht_hq = make_H_tang_hquad(points, faces, normals, deltaHn, H0; gaussorder = 5)
Ht_skip = make_H_tang_skip_singul(points, faces, normals, deltaHn, H0; gaussorder = 5)
psi = PotentialSimple(points, faces, mu, H0; normals = normals)
Ht_erd = HtField(points, faces, psi, normals)
Ht_erd = sqrt.(sum(Ht_erd .* Ht_erd, dims = 1))

Hteor = field_theor(1, 1, 1, mu, H0)
Hn_teor = zeros(size(points,2))
Ht_teor = zeros(size(points,2))
Ht_lang = zeros(size(points,2))
tangs = zeros(size(normals))
for ykey in 1:size(points,2)
    Hn_teor[ykey] = dot(Hteor,normals[:,ykey])
    Ht_teor[ykey] = norm(cross(Hteor,normals[:,ykey]))

    Ht_lang[ykey] = norm(Hteor - Hn_teor[ykey] * normals[:,ykey])
    tangs[:,ykey] = (Hteor - Hn_teor[ykey] * normals[:,ykey]) / norm(Hteor - Hn_teor[ykey] * normals[:,ykey])
end


ENV["MPLBACKEND"]="tkagg"
using PyPlot

pygui(true)
fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")
(x, y, z) = [points[i,:] for i in 1:3]
(vx, vy, vz) = [(tangs .* Ht_teor')[i,:] for i in 1:3]
#ax[:quiver](x,y,z,vx,vy,vz, length=0.6, arrow_length_ratio=0.5, color=:black)
(vx, vy, vz) = [(tangs .* Ht_teor')[i,:] for i in 1:3]
#ax[:quiver](x,y,z,vx,vy,vz, length=0.6, arrow_length_ratio=0.5, color=:blue)
(vx, vy, vz) = [(tangs .* Ht_teor')[i,:] for i in 1:3]
ax[:quiver](x,y,z,vx,vy,vz, length=0.6, arrow_length_ratio=0.5, color=:red)
ax[:scatter](x,y,z, s=2,color="k")
#ax[:quiver](x,y,z,vx,vy,vz, length=50, arrow_length_ratio=0.5)
ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()
