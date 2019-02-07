using LinearAlgebra
using CSV
using JLD2
using ElTopo

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./stabilization.jl")
include("./functions.jl")


#points_csv= CSV.read("/home/lai/Dropbox/dokt/code/matlab/vertices.csv", header=0)
#faces_csv = CSV.read("/home/lai/Dropbox/dokt/code/matlab/triangles.csv", header=0)


points_csv= CSV.read("./meshes/points_sphere.csv", header=0)
faces_csv = CSV.read("/home/lai/Dropbox/dokt/code/matlab/faces_sphere.csv", header=0)

println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')

points = points .* (4.9 * 10^-1)


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



#
# Hn_teor = zeros(size(points,2))
# Ht_teor = zeros(3, size(points,2))
#
# for i in 1:size(points,2)
#     r = norm(points[:,i])
#     r_xy = norm(points[1:2,i])
#     (x,y,z) = [points[j,i] for j in 1:3]
#     Hn_teor[i] = 3 * norm(H0) * points[3,i] / r / (mu+2)
#     Ht_teor[1,i] = - 3 * norm(H0) * r_xy / r * x*z /r/r_xy/(mu+2)
#     Ht_teor[2,i] = - 3 * norm(H0) * r_xy / r * y*z / r/r_xy/(mu+2)
#     Ht_teor[3,i] = 3 * norm(H0) * r_xy / r * r_xy / r / (mu+2)
# end


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


#
# (Htx, Hty, Htz) = Ht[1,:], Ht[2,:], Ht[3,:]
# #velocitiesn = all_normals .* velocitiesn_norm'
# (vnx, vny, vnz) = [velocitiesn[i,:] for i in 1:3]
# (Hnx, Hny, Hnz) = [Hn[i,:] for i in 1:3]
# Hn_teor = normals .* Hn_teor'
# #Ht_teor = normals .* Ht_teor'
# (Hnx_teor, Hny_teor, Hnz_teor) = [Hn_teor[i,:] for i in 1:3]
# (Htx_teor, Hty_teor, Htz_teor) = [Ht_teor[i,:] for i in 1:3]
# temp = normals .* tensorn
# (fx, fy, fz) = [temp[i,:] for i in 1:3]
# (nx, ny, nz) = [normals[i,:] for i in 1:3]
# (Hx, Hy, Hz) = (Htx+Hnx, Hty+Hny, Htz+Hnz)
