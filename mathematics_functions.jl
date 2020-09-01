function gauss_nonsingular(f::Function, r1,r2,r3,gaussorder)
    # requires: using FastGaussQuadrature
    # gaussian integrtation over a triangle with vertices r1,r2,r3
    # f - function of r
    # gaussorder - order of the gaussian integration
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
    # requires: using FastGaussQuadrature
    # gaussian integrtation over a triangle with vertices r1,r2,r3
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
    # requires: using Cubature
    # integral of q(r)/|r-r1|, where q(r) is nonsingular vector function
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

function gauss_flat_pol_vec(f::Function,r1,r2,r3,gaussorder)
    # requires: using FastGaussQuadrature
    # Gaussian integral over a flat triangle with vertices r1,r2,r3
    # integral of f(r), could be 1/r singular
    # singularity on r1

    n_triang = cross(r3-r2,r3-r1) / norm(cross(r3-r2,r3-r1))

    r1_loc, r2_loc, r3_loc = r1-r1, r2-r1, r3-r1
    r2_loc, r3_loc = to_local(r2_loc,n_triang), to_local(r3_loc,n_triang)
    x2, x3 = r2_loc[1], r3_loc[1]
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

    if theta2>theta3 # change of sign if integration limits were reversed
        intval *= -1
    end

    return intval
end

function gauss_curved_pol_vec(f::Function,r1,r2,r3,n1,CDE,gaussorder)
    # requires: using FastGaussQuadrature
    # gaussian integral over a "curved" triangle that is described with a paraboloid
    # with an origin at r1 in the local coordinates of r1 described as:
    # z = C * x^2 + D * x * y + E*y^2
    # z is in the direction of normal at r1: n1
    # integral of f(r), there can be a 1/r singularity at r1
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
