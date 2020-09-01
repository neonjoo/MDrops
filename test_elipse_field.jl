using QuadGK
using CSV
using Plots
using PyPlot

pygui()

include("./functions.jl")
include("./mesh_functions.jl")

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


#points_csv= CSV.read("./meshes/points_ellipse_fewN.csv", header=0)
#faces_csv = CSV.read("./meshes/faces_ellipse_fewN.csv", header=0)

println("Loaded mesh")

# points = convert(Array, points_csv)
# faces = convert(Array, faces_csv)
# points = Array{Float64}(points')
# faces = Array{Int64}(faces')
@load "./meshes/points_critical_hyst_2_21.jld2"
@load "./meshes/faces_critical_hyst_2_21.jld2"

#points, faces = expand_icosamesh(R=1, depth=2)
points = Array{Float64}(points')
faces = Array{Int64}(faces')
normals = Normals(points, faces)

H0 = [0., 0., 1.]
a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])
#println(a, b, c)
#a,b,c = 1, 1, 2.21
mu = 2
mu0 = 4*pi*10e-7
H0 = [0,0,1]
normals = Normals(points, faces)

Hin_teor = field_theor(a,b,c,mu,H0)
H_teor = zeros(3, size(points,2))
Hn_teor = sum(normals .* Hin_teor, dims=1)
for i in 1:size(points,2)
    H_teor[:,i] = Hin_teor .* [1,1,1]
end


psi = PotentialSimple(points, faces, mu, H0; normals = normals)
Ht = HtField(points, faces, psi, normals)
Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = tnormals)
Hn = normals .* Hn_norms'

H_num = Ht+Hn


res = sum((H_num - H_teor).^2, dims=1)


fig = figure(figsize=(7,7))
ax = fig[:gca](projection="3d")
(x, y, z) = [points[i,:] for i in 1:3]
#(vnx, vny, vnz) = [velocities[i,:] for i in 1:3]
ax[:scatter](x,y,z, s=2,color="k")
H_num = H_gauss
#ax[:quiver](x,y,z,H_num[1],H_num[2],H_num[3], length=1, arrow_length_ratio=0.5, color="black")
#ax[:quiver](x,y,z,H_teor[1],H_teor[2],H_teor[3], length=3000, arrow_length_ratio=0.5, color="red")
ax[:set_xlim](-2,2)
ax[:set_ylim](-2,2)
ax[:set_zlim](-2,2)
ax[:set_xlabel]("x axis")
ax[:set_ylabel]("y axis")
ax[:set_zlabel]("z axis")
fig[:show]()
