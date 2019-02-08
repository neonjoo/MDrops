using LinearAlgebra
using CSV
using Makie
#using JLD2
#using ElTopo
#using PyPlot



include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./functions.jl")

points_csv= CSV.read("./meshes/points_sphere.csv", header=0)
faces_csv = CSV.read("./meshes/faces_sphere.csv", header=0)

println("Loaded mesh")

points = convert(Array, points_csv)
faces = convert(Array, faces_csv)
points = Array{Float64}(points')
faces = Array{Int64}(faces')


# H0 = [0, 0, 10]
# eta = 1
# mu = 30
# gamma = 6.9 * 10^-1
# w = 2*pi/50
#
# dt = 0.01
# steps = 50
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
steps = 20



for i in 1:steps

    global points, faces
    global scene
    normals = Normals(points, faces)

    psi = PotentialSimple(points, faces, mu, H0; normals = normals)
    Ht = HtField(points, faces, psi, normals)
    Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
    Hn = normals .* Hn_norms'

    mup = mu
    # magnitudes squared of the normal force
    Hn_2 = sum(Hn.^2, dims=1)
    # magnitudes squared of the tangential force
    Ht_2 = sum(Ht.^2, dims=1)

    tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2

    # make the force normal to surface
    #tensorn = normals .* tensorn

    velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)

    velocitiesn = normals .* velocitiesn_norms'

    points += velocitiesn * dt

    println("first points = $(points[:,1])")

end

scene = Makie.mesh(points', faces',color = :white, shading = false)
Makie.wireframe!(scene[end][1], color = :black, linewidth = 1)
