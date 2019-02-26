datadir="/home/andris/mydatadirst_tuesday_erdmanstabilitytest//"
file = readdir(datadir)[42]
@load "$datadir/$file" data
points, faces, dt = data[1], data[2], data[3]
faces = Array{Int64,2}(faces)
normals = Normals(points, faces)
edges = make_edges(faces)
connectivity = make_connectivity(edges)
println("before = ",minimum(sum(x->x!=0,connectivity,dims=1)))
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)

psi = PotentialSimple(points, faces, mu, H0; normals = normals)
#psi2 = PotentialSimple(points2, faces, mu, H0; normals = normals2)
Ht = HtField(points, faces, psi, normals)
#Ht2 = HtField(points2, faces, psi2, normals2)
Hn_norms = NormalFieldCurrent(points, faces, Ht, mu, H0; normals = normals)
#Hn_norms2 = NormalFieldCurrent(points2, faces, Ht2, mu, H0; normals = normals2)
Hn = normals .* Hn_norms'
#Hn2 = normals2 .* Hn_norms2'

mup = mu
# magnitudes squared of the normal force
Hn_2 = sum(Hn.^2, dims=1)
#Hn2_2 = sum(Hn2.^2, dims=1)
# magnitudes squared of the tangential force
Ht_2 = sum(Ht.^2, dims=1)
#Ht2_2 = sum(Ht2.^2, dims=1)

tensorn = mup*(mup-1)/8/pi * Hn_2 + (mup-1)/8/pi * Ht_2
#tensorn2 = mup*(mup-1)/8/pi * Hn2_2 + (mup-1)/8/pi * Ht2_2

# make the force normal to surface
#tensorn = normals .* tensorn

#velocitiesn_norms = InterfaceSpeedZinchenko(points, faces, tensorn, eta, gamma, normals)
#velocitiesn_norms2 = InterfaceSpeedZinchenko(points2, faces, tensorn2, eta, gamma, normals2)

#velocitiesn = normals .* velocitiesn_norms'
#velocitiesn2 = normals2 .* velocitiesn_norms2'


#velocities = velocitiesn
#velocities2 = make_Vvecs_conjgrad(normals,faces, points, velocitiesn, 1e-6, 120)

velocities = make_magvelocities(points, normals, lambda, Bm, mu, Hn_2, Ht_2)
velocities = sum(velocities .* normals,dims=1) .* normals

zc = SG.Zinchenko2013(points, faces, normals)
SG.stabilise!(velocities,points, faces, normals, zc)

dt = 0.3*minimum(make_min_edges(points,connectivity)./sum(sqrt.(velocities.^2),dims=1))
#dt2 = 0.4*minl2/maxv2
println("dt = $(dt)")
#println("dt2 = $(dt2)")

points += velocities * dt
#points2 += velocities2 * dt
data = [points, faces]
#@save "$datadir/almost_bad2.jld2" data
do_active = false
do_active = flip_edges!(faces, connectivity, points)
println("after = ",minimum(sum(x->x!=0,connectivity,dims=1)))
edges = make_edges(faces)
connectivity = make_connectivity(edges)
println("giga after = ",minimum(sum(x->x!=0,connectivity,dims=1)))

