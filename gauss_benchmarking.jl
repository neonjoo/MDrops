# testing on the cases outlined in
#"Numerical quadrature over smooth, closed surfaces" by J. A. Reeger et al.
#https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2016.0401
using CSV
using FastGaussQuadrature
using LinearAlgebra
using StatsBase
cd("/home/andris/MDrops/")
include("./stabilization.jl")
include("./functions.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")

points_csv= CSV.read("./meshes/points_peanut_a0.25639b0.32049moreN.csv", header=0)
faces_csv = CSV.read("./meshes/faces_peanut_a0.25639b0.32049moreN.csv", header=0)
a = 0.25639
b = 0.32049
meshlambda = a/b
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

function gauss_over_surface(f::Function,points,faces,normals,gaussorder)
    I = 0.
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
        n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]


    end
    return I
end

# testing volume integral
# f = 1/3 * r * n
function gauss_volint(points,faces,normals,gaussorder)
    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end

    function make_f(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        return 1/3 * dot(x, make_n(x,x1,x2,x3,n1,n2,n3))
    end

    I = 0
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
        n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

        I += gauss_nonsingular(x->make_f(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
    end
    return I
end

function trap_volint(points,faces,normals)
    I = 0
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
        n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

        dS = 0.5 * norm(cross(x2-x1,x3-x1))

        f1 = 1/3 * dot(x1,n1)
        f2= 1/3 * dot(x2,n2)
        f3 = 1/3 * dot(x3,n3)

        I += 1/3 * (f1+f2+f3) * dS
    end
    return I
end

function true_volint(a,b)
    return pi/(6a) * (2a*(b^2-2a^2)*sqrt(a^2+b^2) + 3b^4 * asinh(2a*sqrt(a^2+b^2)/b^2))
end

# testing steep gradient function
# f = (2/pi) * atan(100 z)
# analytically: int(f) = 0
function gauss_steepint(points,faces,gaussorder)
    I = 0
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]

        I += gauss_nonsingular(x->2/pi * atan(100 * x[3]), x1,x2,x3,gaussorder)
    end
    return I
end

function trap_steepint(points,faces)
    I = 0
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]

        dS = 0.5 * norm(cross(x2-x1,x3-x1))

        f1 = 2/pi * atan(100 * x1[3])
        f2= 2/pi * atan(100 * x2[3])
        f3 = 2/pi * atan(100 * x3[3])

        I += 1/3 * (f1+f2+f3) * dS
    end
    return I
end

# testing sng function
# f = sign(z)
# analytically: int(f) = 0
function gauss_sgnint(points,faces,gaussorder)
    I = 0
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]

        I += gauss_nonsingular(x->sign(x[3]), x1,x2,x3,gaussorder)
    end
    return I
end

function trap_sgnint(points,faces)
    I = 0
    for i = 1:size(faces,2)
        x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]

        dS = 0.5 * norm(cross(x2-x1,x3-x1))

        f1 = sign(x1[3])
        f2= sign(x2[3])
        f3 = sign(x3[3])

        I += 1/3 * (f1+f2+f3) * dS
    end
    return I
end
