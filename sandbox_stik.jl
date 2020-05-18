using CSV
using FastGaussQuadrature
using LinearAlgebra

#points_csv= CSV.read("./meshes/points_ellipse_manyN.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_ellipse_manyN.csv", header=0)
#println("Loaded mesh")

#points = convert(Array, points_csv)
#faces = convert(Array, faces_csv)
#points = Array{Float64}(points')
#faces = Array{Int64}(faces')

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

# center of mass times S for equilateral triangle with side 2, should be (0,0,0)
println("weak")
println(gauss_weaksingular(r->r*norm(r-[0,2/3*sqrt(3),0]),[0,2/3*sqrt(3),0],[1,-sqrt(3)/3,0],[-1,-sqrt(3)/3,0],5))
# triangle area for equilateral triangle with side 2, should be sqrt(3)=1.732
println(gauss_weaksingular(r->1*norm(r-[0,2/3*sqrt(3),0]),[0,2/3*sqrt(3),0],[1,-sqrt(3)/3,0],[-1,-sqrt(3)/3,0],5))
println("strong")
# center of mass times S for equilateral triangle with side 2, should be (0,0,0)
println(gauss_nonsingular(r->r,[0,2/3*sqrt(3),0],[1,-sqrt(3)/3,0],[-1,-sqrt(3)/3,0],2))
# triangle area for equilateral triangle with side 2, should be sqrt(3)=1.732
println(gauss_nonsingular(r->1,[0,2/3*sqrt(3),0],[1,-sqrt(3)/3,0],[-1,-sqrt(3)/3,0],2))
println()
