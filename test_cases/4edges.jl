
include("../SurfaceGeometry/dt20L/src/SurfaceGeometry.jl")
SG = SurfaceGeometry
#include("./SurfaceGeometry/dt20L/src/StabilisationMethods/stabilisationV2.jl")

#include("./SurfaceGeometry/dt20L/src/Properties.jl")

#include("./SurfaceGeometry/dt20L/src/Utils.jl")

include("../stabilization.jl")
include("../functions.jl")
include("../mesh_functions.jl")
include("../physics_functions.jl")


datadir="./test_cases/"

H0 = [0.,0.,1.]
mu = 30.
Bm = 25.
lambda = 10.


file = readdir(datadir)[2]
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

points0 = copy(points)
faces0 = copy(faces)

@save "$datadir/almost_bad.jld2" data
do_active = false
do_active = flip_edges!(faces, connectivity, points)
println("after = ",minimum(sum(x->x!=0,connectivity,dims=1)))
edges = make_edges(faces)
connectivity = make_connectivity(edges)
println("giga after = ",minimum(sum(x->x!=0,connectivity,dims=1)))


fours = []
for i in 1:size(con1,2)

    clean = filter(x -> x!=0, con1[:,i])
    if size(clean,1) == 4
        push!(fours, i)
    end
end





function flip_edges2!(faces, connectivity, vertices)
    # flip edges to improve mesh, updates "faces" and "connectivity"

    println("flipping edges")
    maxx = size(vertices, 2)
    global continue_flip = true
    global flipped_any = false

    while continue_flip
        global continue_flip
        continue_flip = false
        #println("started flip loop")

        for i in fours
            println("------------------------ checking $i in fours")
            # num of i-th vertex neighbors
            i_num = length(filter(x -> x>0, connectivity[:,i]))
            println("^ i_num: $i_num")
            if i_num <= 5
                continue
            end

            for j in filter(x -> x>0, connectivity[:,i])
                println("-- checking $i-$j")
                # num of j-th vertex neighbors
                j_num = length(filter(x -> x>0, connectivity[:,j]))
                println("$i-$j: j_num: $j_num")
                if j_num <= 5
                    continue
                end

                xi, xj = vertices[:,i], vertices[:,j]

                common = intersect(connectivity[:,i], connectivity[:,j])
                common = filter(x -> x>0, common)

                if length(common) != 2
                    continue
                end

                k, m = common[1], common[2]
                xk, xm = vertices[:,k], vertices[:,m]

                kc = find_circumcenter(xi, xj, xk)
                mc = find_circumcenter(xi, xj, xm)

                # a threshold value for flipping (Zinchenko2013)
                d = norm(dot(xk-kc, xm-xk)) + norm(dot(xm-mc, xm-xk))

                if norm(xk - xm)^2 < d
                    println("--- flippening $i--$j to $k--$m")
                    flip_connectivity2!(faces, connectivity, i, j, k, m)
                    continue_flip = true
                    break
                end

            end # end j for
            if continue_flip
                flipped_any = true
                break
            end
        end # end i for

    end # end while

    # returns true if any edge was flipped. If so, active stabilization is to be applied
    println("--------- Flipped any?  $flipped_any ---- ")
    return flipped_any

end # end function


function flip_connectivity2!(faces, connectivity, i, j, k, m)
    # adjusts faces & connectivity to the i--j  ->  k--m edge flip

    found_one, found_two = false, false
    for s in 1:size(faces,2)
        # finds first(and only) column where all 3 indices appear and adjusts indices

        if !found_one
            if length(intersect([i,j,k], faces[:,s])) == 3
                row = findfirst(x-> x==j, faces[:,s])
                faces[row, s] = m
                found_one = true
            end
        end

        if !found_two
            if length(intersect([i,m,j], faces[:,s])) == 3
                row = findfirst(x-> x==i, faces[:,s])
                faces[row, s] = k
                found_two = true
            end
        end

        if found_one && found_two
            break
        end

    end # end s for


    # row of j in i-th column etc.
    row_j_in_i = findfirst(x-> x==j, connectivity[:,i])
    row_i_in_j = findfirst(x-> x==i, connectivity[:,j])
    row_k_in_m = findfirst(x-> x==0, connectivity[:,m])
    row_m_in_k = findfirst(x-> x==0, connectivity[:,k])

    # cut the i--j edge in "connectivity" by sliding column values up 1 row, placing 0 in the end
    connectivity[row_i_in_j:end, j] = [connectivity[row_i_in_j+1:end, j]; 0]
    connectivity[row_j_in_i:end, i] = [connectivity[row_j_in_i+1:end, i]; 0]

    # adds a zero row if either of new vertices k or m are connected to some previous maximum connectivity (aka valence)
    if row_k_in_m == nothing || row_m_in_k == nothing
        println("padded row of 0s on bottom of \"connectivity\"")

        connectivity = vcat(connectivity, zeros(Int64,1, size(connectivity,2)))

        if row_k_in_m == nothing
            connectivity[end, m] = k
        else
            connectivity[row_k_in_m, m] = k
        end

        if row_m_in_k == nothing
            connectivity[end, k] = m
        else
            connectivity[row_m_in_k, k] = m
        end

    else
        connectivity[row_k_in_m, m] = k
        connectivity[row_m_in_k, k] = m
    end # end if


end # end function
