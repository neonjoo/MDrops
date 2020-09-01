using Pkg
pkg"activate ."
pkg"resolve"

using Plots
using StatsBase
using Optim
#p = Plots
#p.pyplot()
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./sandbox_lang.jl")

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

function gauss_curved(f::Function,r1,r2,r3,n1,CDE,gaussorder)
    # integral of f(r)
    r2 = project_on_given_paraboloid(CDE, n1, r2, r1)
    r3 = project_on_given_paraboloid(CDE, n1, r3, r1)
    r1_loc, r2_loc, r3_loc = r1-r1, r2-r1, r3-r1
    r2_loc, r3_loc = to_local(r2_loc,n1), to_local(r3_loc,n1)

    x2, x3 = r2_loc[1], r3_loc[1] + 1e-15
    y2, y3 = r2_loc[2], r3_loc[2]

    # angles are from -pi to pi, masured from the x axis
    theta2, theta3 = atan(y2,x2), atan(y3,x3)


    # effectively change integration limits so that theta sweeps only the triangle
    if (theta3 - theta2) > pi
        theta2 += 2pi
    end

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

    function make_r(rho,theta, n1, r1, CDE)
        x = rho * cos(theta)
        y = rho * sin(theta)
        z = CDE[1]*x^2 + CDE[2]*x*y + CDE[3]*y^2

        r_loc = [x, y, z]

        return to_global(r_loc,n1) + r1
    end

    intval = [0.,0.,0.]

    for i = 1:length(v)
        theta = make_theta(v[i], theta2, theta3)
        rho_m = make_rho_m(theta, x2, x3, y2, y3)
        for k = 1:length(u)
            rho = make_rho(u[k], rho_m)
            x, y = rho*cos(theta), rho*sin(theta)
            dA = sqrt( (2*CDE[1]*x + CDE[2]*y)^2 + (CDE[2]*x + 2*CDE[3]*y)^2 + 1)
            intval += wv[i]*(theta3 - theta2)/2 *
                    wu[k]*f( make_r(rho,theta, n1, r1, CDE) ) * rho * rho_m/2 * dA
        end
    end
    if theta2 > theta3 # change of sign if integration limits are reversed
        return -intval
    end

    return intval
end

function make_H_field(points, faces, normals, mu, H0, gaussorder)

    N = size(points, 2)
    Ht = zeros(3, N) # 3xN array
    deltaS = make_dS(points,faces)

    edges = make_edges(faces)
    connectivity = make_connectivity(edges)
    normals, CDE = make_normals_spline(points, connectivity, edges, normals)
    dHn = make_deltaH_normal(points, faces, mu, H0)

    for ykey in 1:N
        for i in 1:size(faces, 2)
            function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                n = B * zeta_xi_eta
                return n/norm(n)
            end

            function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta normal field linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [dHn1 dHn2 dHn3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                deltaHn = dot(B, zeta_xi_eta)
                return deltaHn
            end

            function make_Ht(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = y - x
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                temp = dHnx * ny - dHn[ykey] * nx
                Ht_func = cross(temp, r) / norm(r)^3

                return Ht_func/(4pi)
            end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1, dHn1 = normals[:,faces[1,i]], dHn[faces[1,i]]
                n2, dHn2 = normals[:,faces[2,i]], dHn[faces[2,i]]
                n3, dHn3 = normals[:,faces[3,i]], dHn[faces[3,i]]

                Ht[:, ykey] += gauss_nonsingular(x->make_Ht(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,gaussorder)
            else # if is singular triangle
                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])
                ind2, ind3 = (singul_ind) % 3 + 1, (singul_ind + 1) % 3 + 1

                x1, x2, x3 = points[:,faces[singul_ind,i]], points[:,faces[ind2,i]], points[:,faces[ind3,i]]

                n1, dHn1 = normals[:,faces[singul_ind,i]], dHn[faces[singul_ind,i]]
                n2, dHn2 = normals[:,faces[ind2,i]], dHn[faces[ind2,i]]
                n3, dHn3 = normals[:,faces[ind3,i]], dHn[faces[ind3,i]]
                Ht[:,ykey] += gauss_curved(x->make_Ht(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,n1,CDE[:,ykey],gaussorder)

            end # end if singular

        end # end i for

        Ht[:, ykey] += cross(normals[:, ykey], H0)
    end # end ykey for

    Ht_norms2 = sum(Ht.^2, dims=1)
    Ht_norms = sqrt.(Ht_norms2)

    return Ht_norms
end


function make_H_field_curved(points, faces, normals, mu, H0, gaussorder)

    N = size(points, 2)
    Ht = zeros(3, N) # 3xN array
    deltaS = make_dS(points,faces)

    edges = make_edges(faces)
    connectivity = make_connectivity(edges)
    normals, CDE = make_normals_spline(points, connectivity, edges, normals)
    dHn = make_deltaH_normal(points, faces, mu, H0)

    for ykey in 1:N
        for i in 1:size(faces, 2)
            function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [n1 n2 n3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                n = B * zeta_xi_eta
                return n/norm(n)
            end

            function make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3) # delta normal field linear interpolation
                A = [x1 x2 x3] # matrix of vertex radiusvecotrs
                B = [dHn1 dHn2 dHn3] # matrix of vertex normals

                zeta_xi_eta = A \ x # find local triangle parameters

                deltaHn = dot(B, zeta_xi_eta)
                return deltaHn
            end

            function make_Ht(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3; y=points[:,ykey],ny=normals[:,ykey])
                nx = make_n(x,x1,x2,x3,n1,n2,n3)
                r = y - x
                dHnx = make_dHn(x,x1,x2,x3,dHn1,dHn2,dHn3)
                temp = dHnx * ny - dHn[ykey] * nx
                Ht_func = cross(temp, r) / norm(r)^3

                return Ht_func/(4pi)
            end

            if !(ykey in faces[:,i]) # if not singular triangle
                x1 = points[:,faces[1,i]]
                x2 = points[:,faces[2,i]]
                x3 = points[:,faces[3,i]]

                n1, dHn1 = normals[:,faces[1,i]], dHn[faces[1,i]]
                n2, dHn2 = normals[:,faces[2,i]], dHn[faces[2,i]]
                n3, dHn3 = normals[:,faces[3,i]], dHn[faces[3,i]]

                Ht[:, ykey] += gauss_curved(x->make_Ht(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,n1,CDE[:,faces[1,i]],gaussorder)
            else # if is singular triangle
                # arrange labels so that singularity is on x1
                # (singul_ind + n - 1) % 3 + 1 shifts index by n
                singul_ind = findfirst(singul_ind->singul_ind==ykey,faces[:,i])
                ind2, ind3 = (singul_ind) % 3 + 1, (singul_ind + 1) % 3 + 1

                x1, x2, x3 = points[:,faces[singul_ind,i]], points[:,faces[ind2,i]], points[:,faces[ind3,i]]

                n1, dHn1 = normals[:,faces[singul_ind,i]], dHn[faces[singul_ind,i]]
                n2, dHn2 = normals[:,faces[ind2,i]], dHn[faces[ind2,i]]
                n3, dHn3 = normals[:,faces[ind3,i]], dHn[faces[ind3,i]]
                Ht[:,ykey] += gauss_curved(x->make_Ht(x,x1,x2,x3,n1,n2,n3,dHn1,dHn2,dHn3), x1,x2,x3,n1,CDE[:,ykey],gaussorder)

            end # end if singular

        end # end i for

        Ht[:, ykey] += cross(normals[:, ykey], H0)
    end # end ykey for

    Ht_norms2 = sum(Ht.^2, dims=1)
    Ht_norms = sqrt.(Ht_norms2)

    return Ht_norms
end


@load "./meshes/points_critical_hyst_2_21.jld2"
@load "./meshes/faces_critical_hyst_2_21.jld2"

#points, faces = expand_icosamesh(R=1, depth=2)
points = Array{Float64}(points')
faces = Array{Int64}(faces')
normals = Normals(points, faces)

H0 = [0., 0., 1.]
a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])

#points, faces = expand_icosamesh(R=1, depth=2)
#points = Array{Float64}(points)
#faces = Array{Int64}(faces)
#normals = Normals(points, faces)

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals, CDE = make_normals_spline(points, connectivity, edges, normals)

H0 = [0., 0., 1.]

ht_elips5 = make_H_field(points, faces, normals, 2, H0, 5)
ht_celips3_2node = make_H_field_curved(points, faces, normals, 2, H0, 5)
hn_elips3 = make_deltaH_normal(points, faces, 2, H0)
