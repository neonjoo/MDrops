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

#points_csv= CSV.read("/home/lai/Dropbox/dokt/code/matlab/vertices.csv", header=0)
#faces_csv = CSV.read("/home/lai/Dropbox/dokt/code/matlab/triangles.csv", header=0)

points_csv= CSV.read("/home/lai/Dropbox/dokt/code/matlab/points_sphere.csv", header=0)
faces_csv = CSV.read("/home/lai/Dropbox/dokt/code/matlab/faces_sphere.csv", header=0)

println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')

points = points .* (4.9 * 10^-1)

# all_normals = zeros(Float64, 3, length(points[1,:]))
# all_normals_theor = all_theor(points, coef)

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

#psi_noreg = PotentialSimple(points, faces, mu, H0; normals = all_normals, regularize=false)

# psi_theor = zeros(Float64, size(points,2))
#
# for i in 1:size(points,2)
#     psi_theor[i] += dot(H0, points[:,i])
# end


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



H0 = [0, 0, 10]


#H0 = 33 .* [0, 0, 1]
#H0 = [0,0,0]
mu = 30
eta = 1
gamma = 6.9 * 10^-1

continue_sim = false
last_step = 0

w = 2*pi/50
t = 0
dt = 0.01
steps = 10

datadir="/home/lai/Dropbox/dokt/code/data/stikuts3"
# typical triangle side length
scale = 0.1 * 4.9 * 10^-1

elparameters(scale) = Elparameters( # comments are some values that work fine
 m_use_fraction = false,                     # false
 m_min_edge_length = 0.7*scale,             # 0.7 * scale
 m_max_edge_length = 1.5*scale,               # 1.5 * scale
 m_max_volume_change = 0.1*scale^3,         # 0.1 * scale^3
 m_min_curvature_multiplier = 1,             # 1
 m_max_curvature_multiplier = 1,            # 1
 m_merge_proximity_epsilon = 0.5*scale,     # 0.5 * scale
 m_proximity_epsilon = 0.00001,             # 0.00001
 m_perform_improvement = true,              # true
 m_collision_safety = false,                 # false
 m_min_triangle_angle = 15,                 # 15
 m_max_triangle_angle = 120,                # 120
 m_allow_vertex_movement = false,           # false   ### This is where is a bug
 m_use_curvature_when_collapsing = false,    # false
 m_use_curvature_when_splitting = false,    # false
 m_dt = 1                                   # 1
)

par = elparameters(scale)


trajectory = zeros(steps + 1, 3, size(points,2))
trajectory[1,:,:] = points

global points2 = copy(points)
global faces2 = copy(faces)

global normals, Ht, Hn, velocitiesn, velocitiesn_norms


if !isdir("$datadir")
    mkdir("$datadir")

    open("$datadir/_params.txt", "w") do file
        write(file, "H0=$H0\nmu=$mu\neta=$eta\ngamma=$gamma\nsteps=$steps")
    end

    println("Created new dir: $datadir")
end



if continue_sim
    last_file = readdir(datadir)[end]
    println(last_file)

    # fix hardcoding
    global data
    @load "$datadir/data02000.jld2" data
    global last_step = 0
end

#@load "/home/lai/Dropbox/dokt/code/data/elongation2/data14000.jld2" data

#global points2, faces2, H0 = data[1], data[2], data[3]

expand = true
reset_vmax = true

for i in 1:steps
    global points2, faces2, t, w, H0
    global points = copy(points2)
    global faces = copy(faces2)
    global expand, reset_vmax

    #normals = zeros(3, size(points,2))
    global normals = Normals(points, faces)

    global psi = PotentialSimple(points, faces, mu, H0; normals = normals)
    global Ht = HtField(points, faces, psi, normals)
    global Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
    global Hn = normals .* Hn_norms'

    mup = mu
    # magnitudes squared of the normal force
    Hn_2 = sum(Hn.^2, dims=1)
    # magnitudes squared of the tangential force
    Ht_2 = sum(Ht.^2, dims=1)


    global tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2

    # make the force normal to surface
    #tensorn = normals .* tensorn

    global velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)

    global velocitiesn = normals .* velocitiesn_norms'



    velocitiesn = make_Vvecs_conjgrad(normals,faces, points, velocitiesn, 1e-6, 120);

    velocitiesn = velocitiesn'

    dt = 0.1 * scale / max(sqrt(sum(Vvecs.*Vvecs,2)))

    points2 = points + velocitiesn * dt
    # trajectory[i + 1,:,:] = points2

    # ElTopo magic
    #actualdt,points2,faces2 = improvemeshcol(points,faces,points2,par)

    if i % 1 == 0
        data = [points2, faces2, (H0)]
        println("Finished step $(last_step + i)")
        @save "$datadir/data$(lpad(i + last_step,5,"0")).jld2" data

    end


    if reset_vmax
        println("Resetting v0max")
        global v0max = maximum(abs.(velocitiesn))
        reset_vmax = false
    end

    vi = maximum(abs.(velocitiesn))

    println("vi = $vi, v0max = $v0max, vi/v0max = $(vi/v0max)")

    if vi/v0max < 0.05
        println("Reached eqlb with vi=$vi")
        global expand, reset_vmax
        reset_vmax = true

        if expand
            println("expandasdas")
            H0 += [0,0,0.5]
        else
            H0 -= [0,0,0.05]
        end
        if H0[3] >= 3.5
            global expand = false
        end
    end


    t += dt
    #H0 = 33 .* [sin(w*t), 0, cos(w*t)]
end

#println(H0, mu, st.mean(velocitiesn_norms))

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

Hn_teor = zeros(size(points,2))
Ht_teor = zeros(3, size(points,2))

for i in 1:size(points,2)
    r = norm(points[:,i])
    r_xy = norm(points[1:2,i])
    (x,y,z) = [points[j,i] for j in 1:3]
    Hn_teor[i] = 3 * norm(H0) * points[3,i] / r / (mu+2)
    Ht_teor[1,i] = - 3 * norm(H0) * r_xy / r * x*z /r/r_xy/(mu+2)
    Ht_teor[2,i] = - 3 * norm(H0) * r_xy / r * y*z / r/r_xy/(mu+2)
    Ht_teor[3,i] = 3 * norm(H0) * r_xy / r * r_xy / r / (mu+2)
end



(Htx, Hty, Htz) = Ht[1,:], Ht[2,:], Ht[3,:]
#velocitiesn = all_normals .* velocitiesn_norm'
(vnx, vny, vnz) = [velocitiesn[i,:] for i in 1:3]
(Hnx, Hny, Hnz) = [Hn[i,:] for i in 1:3]
Hn_teor = normals .* Hn_teor'
#Ht_teor = normals .* Ht_teor'
(Hnx_teor, Hny_teor, Hnz_teor) = [Hn_teor[i,:] for i in 1:3]
(Htx_teor, Hty_teor, Htz_teor) = [Ht_teor[i,:] for i in 1:3]
temp = normals .* tensorn
(fx, fy, fz) = [temp[i,:] for i in 1:3]
(nx, ny, nz) = [normals[i,:] for i in 1:3]
(Hx, Hy, Hz) = (Htx+Hnx, Hty+Hny, Htz+Hnz)

using PyPlot
pygui()

fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")

(x, y, z) = [points[i,:] for i in 1:3]
ax[:scatter](x,y,z, s=2,color="k")
#(xp, yp, zp) = points[:,418]
#ax[:scatter](xp, yp, zp, s=10,color="r")

#ax[:set_title]("Normals")
#ax[:quiver](x,y,z,nx,ny,nz, length=0.2, arrow_length_ratio=0.5)

#ax[:quiver](x,y,z,fx,fy,fz, length=0.2, arrow_length_ratio=0.5)
#ax[:set_title]("Normal forces")

ax[:quiver](x,y,z,vnx,vny,vnz, length=5, arrow_length_ratio=0.5)
#ax[:set_title]("Normal velocities")

#ax[:quiver](x,y,z,Hnx,Hny,Hnz, length=0.3, arrow_length_ratio=0.5, color="red")
#ax[:quiver](x,y,z,Hnx_teor,Hny_teor,Hnz_teor, length=0.3, arrow_length_ratio=0.5,color="blue")
#ax[:set_title]("Normal H component, external H // z")

#ax[:quiver](x,y,z,Htx,Hty,Htz, length=0.3, arrow_length_ratio=0.5, color="blue")
#ax[:quiver](x,y,z,Htx_teor,Hty_teor,Htz_teor, length=0.3, arrow_length_ratio=0.5,color="red")
#ax[:set_title]("Tangential H")

#ax[:quiver](x,y,z,Hx,Hy,Hz, length=0.3)
ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()



function make_trianpath(vertices,triangles,i)
    trianpath = [
        vertices[triangles[i,1],:]
        vertices[triangles[i,2],:]
        vertices[triangles[i,3],:]
        vertices[triangles[i,1],:]
    ]

    return transpose(reshape(trianpath,3,4))
end

using PyPlot

using Plots

function mesh_plot(vertices,triangles)
    #fig = figure(figsize(7,7))
    #closeall()

    #fig = figure(figsize=(7,7))
    #ax = fig[:gca](projection="3d")

    #(x, y, z) = [points2[i,:] for i in 1:3]
    #ax[:scatter](x,y,z, s=2,color="k")
    plotly()
    vertices = vertices'
    triangles = triangles'
    points = make_trianpath(vertices,triangles,1)
    p = plot(points[:,1],points[:,2],points[:,3],c="black")
    for i in 1:20#size(triangles)[1]
        points = make_trianpath(vertices,triangles,i)
        plot!(p,points[:,1],points[:,2],points[:,3],c="black")
    end
    #plot!(xlims=(-2,2),ylims=(-2,2),zlims=(-2,2),legend=false)
    println("Done. Showing..")
    #gui()
    show(p)
end
