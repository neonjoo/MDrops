cd("/home/andris/MDrops/")

using Pkg

pkg"activate ."
pkg"resolve"

using JLD2
using StatsBase
using LinearAlgebra
using FastGaussQuadrature
using Optim
using Plots

include("./SurfaceGeometry/dt20L/src/Iterators.jl")
include("./mesh_functions.jl")
include("./physics_functions.jl")
include("./mathematics_functions.jl")
## making mesh

datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_adaptiveN_adaptive_dt_old_surface_stabil_flip2_splitfrom13/"

files = readdir(datadir)
N = 80
file = files[2+N]
println(file)
Ndata = size(files,1)-3
# file = files[4+N]
# file = files[525]
@load "$datadir/$file" data

points, faces = data[1], data[2]
faces = Array{Int64,2}(faces)
mean_x, mean_y, mean_z = (StatsBase.mean(points[1,:]),
                StatsBase.mean(points[2,:]),
                StatsBase.mean(points[3,:]))
points = points .- [mean_x, mean_y, mean_z]


a,b,c = maximum(points[1,:]), maximum(points[2,:]), maximum(points[3,:])

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)
(normals, CDE) = make_normals_spline(points, connectivity, edges, normals)


k1, k2 = make_pc(CDE)
H = sqrt.(k1.^2 + k2.^2)
deltaS = make_dS(points,faces)
Hface = zeros(size(faces,2))
deltaSface = zeros(size(faces,2))
zface = zeros(size(faces,2))
maxLface = zeros(size(faces,2))
for face_ind in 1:size(faces,2)

    for i = 1:3
        v1 = points[:,faces[i,face_ind]]
        v2 = points[:,faces[i% 3 + 1,face_ind] ]
        maxLface[face_ind] = max(maxLface[face_ind], norm(v1-v2))
    end

    Hface[face_ind] = sum(H[faces[:,face_ind]]) / 3
    deltaSface[face_ind] = sum(deltaS[faces[:,face_ind]]) / 3
    zface[face_ind] = sum(points[3,faces[:,face_ind]]) / 3
end
deltaS = make_dS(points,faces)
critS = sqrt.(deltaSface) .* Hface
critL = maxLface .* Hface
## add points to places of high curvature
Plots.scatter(zface,critS)
Plots.scatter!(zface,critL)
#
# function make_neighbor_faces(faces)
#     neighbor_faces = Array{Int64}(undef, 3, 0)
#     for i = 1:size(faces,2)
#
#         vind_1, vind_2, vind_3 = faces[:,i] # vertex indices
#         edge1, edge2, edge3 = sort([vind_1,vind_2]), sort([vind_1,vind_3]), sort([vind_2,vind_3]) #edges
#         edges_i = [edge1,edge2,edge3]
#         # find neighbors for triangle i
#         this_neighbors = [0,0,0]
#         this_neighbors_ind = 1
#         for j = 1:size(faces,2)
#             if i != j
#                 vjind_1, vjind_2, vjind_3 = faces[:,j] # vertex indices
#                 edgej1, edgej2, edgej3 = sort([vjind_1,vjind_2]), sort([vjind_1,vjind_3]), sort([vjind_2,vjind_3]) #edges
#
#                 # if any of the edges of j are in i
#                 if (edgej1 in edges_i) || (edgej2 in edges_i) || (edgej3 in edges_i)
#                     this_neighbors[this_neighbors_ind] = j
#                     this_neighbors_ind += 1
#                     if this_neighbors_ind == 5
#                         #this_neighbors_ind = 1
#                         println("hmm, somtheing wrong - 4 neighbors")
#                     end
#                 end
#             end
#         end
#
#         neighbor_faces = cat(neighbor_faces,this_neighbors,dims=2)
#     end
#     return neighbor_faces
# end
#
# function mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = 0.55)
#     k1, k2 = make_pc(CDE) # principal curvatures on vertices
#     H = sqrt.(k1.^2 + k2.^2)
#
#     ## mark faces that are too curved and need to be split
#     marked_faces = falses(size(faces,2))
#
#     for i = 1:size(faces,2)
#         v1, v2, v3 = points[:,faces[1,i]], points[:,faces[2,i]], points[:,faces[3,i]] # triangle vertices
#         d1, d2, d3 = norm(v1-v2), norm(v1-v3), norm(v2-v3) # edge lengths
#
#         maxd = max(d1,d2,d3) # longest edge
#         Hface = sum(H[faces[:,i]]) / 3 # average curvature H of vertices
#
#         crit = Hface * maxd
#         if crit > cutoff_crit
#             marked_faces[i] = true
#         end
#     end
#
#     # mark also all the faces that have at least two nearby marked faces
#     neighbors_marked = false
#     while !neighbors_marked
#         marked_faces_old = copy(marked_faces)
#         for i = 1:size(faces,2)
#             if marked_faces[i] == false
#
#                 nearby_marked_faces = 0
#                 for j in neighbor_faces[:,i]
#                     if marked_faces[j] == true
#                         nearby_marked_faces += 1
#                     end
#                 end
#                 if nearby_marked_faces > 1
#                     marked_faces[i] == true
#                 end
#
#             end
#         end
#
#         if marked_faces_old == marked_faces
#             neighbors_marked = true
#         end
#     end # end while
#
#
#     return marked_faces
# end
#
# neighbor_faces = make_neighbor_faces(faces)
#
#
#
# function add_points(points, faces,normals, edges, CDE; cutoff_crit = 0.55)
#     # add points in places where the curvature is too great
#     println("started adding points")
#     neighbor_faces = make_neighbor_faces(faces)
#     marked_faces = mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = cutoff_crit)
#     println("Number of marked triangles: ",sum(marked_faces))
#
#     faces_vec = [faces[:,i] for i in 1:size(faces,2)] # makes it easier to delete or add a face
#     new_points = copy(points)
#     checked_sides = Array{Int64}[] # create empty array
#     is = Array{Int64}(undef,(1,3))
#     js = Array{Int64}(undef,(1,3))
#     ks = Array{Int64}(undef,(1,3))
#     original_vertices_length = size(points,2)
#
#     for i in 1:size(faces,2)
#     if marked_faces[i] == true # loop throough face to be split
#         #println("splitting marked triangle #$i")
#         triangle = faces[:,i]
#
#         is[1] = triangle[1]
#         is[2] = triangle[2]
#         is[3] = triangle[3]
#
#         sides = [
#             sort([triangle[1], triangle[2]]),
#             sort([triangle[2], triangle[3]]),
#             sort([triangle[3], triangle[1]])
#         ]
#
#         for j = 1:3
#             if sides[j] in checked_sides
#                 js[j] = findfirst(x -> x==sides[j],checked_sides) + original_vertices_length
#             else
#                 new_vertex = 0.5*( new_points[:,sides[j][1]] + new_points[:,sides[j][2]] )
#
#                 proj1 = project_on_given_paraboloid(CDE[:,sides[j][1]],normals[:,sides[j][1]], new_vertex, new_points[:,sides[j][1]])
#                 proj2 = project_on_given_paraboloid(CDE[:,sides[j][2]],normals[:,sides[j][2]], new_vertex, new_points[:,sides[j][2]])
#
#                 new_vertex = (proj1 + proj2)/2
#
#                 push!(checked_sides, sides[j]) # add side to checked sides
#                 new_points = hcat(new_points,new_vertex) # add to new vertices
#                 js[j] = size(checked_sides,1) + original_vertices_length
#             end
#         end
#
#         # remove marked triangle from triangle array
#         filter!(x->x!=triangle,faces_vec)
#
#         # add triangles that the marked triangle was split into
#         push!(faces_vec, [is[1],js[1],js[3]])
#         push!(faces_vec, [js[1],is[2],js[2]])
#         push!(faces_vec, [js[2],is[3],js[3]])
#         push!(faces_vec, [js[1],js[2],js[3]])
#
#         # loop through neighbor triangles and check if they are also marked
#         # if not, split them in two
#
#         for k in neighbor_faces[:,i]
#         if marked_faces[k] == false
#             # identify the point across the marked triangle
#             neigh_side = [] # if neigh_point in triangle
#             neigh_triangle = faces[:,k]
#             for (neigh_ind, neigh_point) in enumerate(faces[:,k])
#                 if neigh_point in triangle
#                     push!(neigh_side,neigh_point)
#                 else
#                     #ks = neigh_point
#                     # shift the indices in the neighboring triangle so that
#                     # the point across the marked triangle is first
#                     neigh_triangle = circshift(faces[:,k],4-neigh_ind)
#                 end
#             end
#             sort!(neigh_side)
#             # point split in the middle of shared vertex
#             jjs = findfirst(x -> x==neigh_side,checked_sides) + original_vertices_length
#             # rename points in this triangle for ease of writing
#             ks = neigh_triangle
#
#
#             # remove adjacent triangle from triangle array
#             filter!(x->x!=faces[:,k],faces_vec)
#             # split it in two
#             push!(faces_vec, [ks[2],jjs,ks[1]])
#             push!(faces_vec, [jjs,ks[3],ks[1]])
#
#         end
#         end # k neighbor triangle loop
#
#     end
#     end # i triangle loop
#
#     new_faces = []
#     for (ind,face) in enumerate(faces_vec)
#         if ind == 1
#             new_faces = face
#         else
#             new_faces = hcat(new_faces,face)
#         end
#     end
#
#     return new_points, new_faces
# end
#
# new_points, new_faces = add_points(points, faces,normals, edges, CDE; cutoff_crit = 0.55)
#
# new_normals = Normals(new_points,new_faces)


function make_deep_normals_spline(points, connectivity, edges, normals0;
                    Cs=1.0, eps_inner=1e-5, eps_outer=1e-4,
                    max_iters_inner=1000, max_iters_outer=100)
    # this function differs from make_normals_spline() by the fact that it takes
    # the vertices that neighbor the neighboring vertices in the paraboloid fitting
    #returns improved normals and CDE - parameters of locally fitted paraboloid
    #z = C*x^2+D*x*y+E*y^2
    #Cs is a coupling parameter. Zinchencko(2000) sets it to 1
    #eps inner and outer are convergence criteria.
    #eps outer >> eps_inner

    # takes almost 1minute to enter this function from main_lan.jl call

    println("fitting local paraboloids, outer costs: ")

    CDE = zeros(Float64, size(points)) # 3xN array
    gradPhi = [0.,0.,0.]
    normals = copy(normals0)
    #outer iterations

    for m = 1:max_iters_outer
        #println(m)
        normalsp = copy(normals) #previous
        for i = 1:size(points,2)
            #inner iterations
            for k = 1:max_iters_inner
                #println("inner iter = ",k)
                # get edges close to vertex
                # ev = 3xn matrix containing edge vectors from i-th to all adjecent points + a level deeper
                # n = v^2-2v
                # ev is padded with 0-os
                v = size(connectivity,1) ## maximum valence
                ev = zeros(Float64,3,v^2-2v)
                edge_normal_sum = zeros(Float64,3,v^2-2v)
                checked_points = zeros(Int64,v^2-2v)

                ind = 1
                for j in connectivity[:,i]
                    # iterate through close points j
                    if j == 0
                        break
                    end

                    if !(j in checked_points)
                        ev[:,ind] = points[:,j] - points[:,i]
                        edge_normal_sum[:,ind] = normals[:,j] + normals[:,i]
                        checked_points[ind] = j

                        ev[:,ind] = to_local(ev[:,ind],normals[:,i])
                        edge_normal_sum[:,ind] = to_local(edge_normal_sum[:,ind],normals[:,i])
                        ind += 1
                    end

                    for k in connectivity[:,j]
                        # iterate through next level of close points j
                        if k == 0
                            break
                        end

                        if !(k in checked_points)
                            #println(checked_points)
                            ev[:,ind] = points[:,k] - points[:,i]
                            edge_normal_sum[:,ind] = normals[:,k] + normals[:,i]
                            checked_points[ind] = k

                            ev[:,ind] = to_local(ev[:,ind],normals[:,i])
                            edge_normal_sum[:,ind] = to_local(edge_normal_sum[:,ind],normals[:,i])
                            ind += 1
                        end

                    end
                end
                # fit parabola and get CDE coefficients

                #M  3x3
                #b  3x1
                # make matrix M and vector b to be solved for least squares problem
                # println(ev[1,:])
                # println(ev[2,:])
                # println(ev[3,:])
                # println(-2*ev[3,:]./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                #     ev[1,:].^2)
                # readline(stdin)

                b = [sum(filter(x -> !isnan(x),
                 -2*ev[3,:]./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^2
                 ));
                 sum(filter(x -> !isnan(x),
                 -2*ev[3,:]./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].*ev[2,:]
                 ));
                 sum(filter(x -> !isnan(x),
                 -2*ev[3,:]./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[2,:].^2
                 ))
                ]

                M = [sum(filter(x -> !isnan(x),
                    2 ./ (ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^4
                     )) sum(filter(x -> !isnan(x),
                    2 ./ (ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^3 .* ev[2,:]
                     )) sum(filter(x -> !isnan(x),
                    2 ./ (ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^2 .* ev[2,:].^2
                     ));


                     sum(filter(x -> !isnan(x),
                    2 ./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^3 .* ev[2,:]
                     )) sum(filter(x -> !isnan(x),
                    2 ./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^2 .* ev[2,:].^2
                     )) sum(filter(x -> !isnan(x),
                    2 ./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:] .* ev[2,:].^3
                     ));


                      sum(filter(x -> !isnan(x),
                    2 ./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:].^2 .* ev[2,:].^2
                     )) sum(filter(x -> !isnan(x),
                    2 ./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[1,:] .* ev[2,:].^3
                     )) sum(filter(x -> !isnan(x),
                    2 ./(ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                     ev[2,:].^4
                     ))
                ]

                CDE[:,i] = M\(-b)

                # get gradPi
                gradPhi[1] = sum(filter(x -> !isnan(x),
                            -2*(CDE[1,i]*ev[1,:].^2 + CDE[2,i]*ev[1,:].*ev[2,:] +
                            CDE[3,i]*ev[2,:].^2 - ev[3,:]) ./
                            (ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                            (2*CDE[1,i]*ev[1,:].*ev[3,:] +
                            CDE[2,i]*ev[2,:].*ev[3,:] + ev[1,:])
                            ))

                gradPhi[2] = sum(filter(x -> !isnan(x),
                            -2*(CDE[1,i]*ev[1,:].^2 + CDE[2,i]*ev[1,:].*ev[2,:] +
                            CDE[3,i]*ev[2,:].^2 - ev[3,:]) ./
                            (ev[1,:].^2 + ev[2,:].^2 + ev[3,:].^2) .*
                            (CDE[2,i]*ev[1,:].*ev[3,:] +
                            2*CDE[3,i]*ev[2,:].*ev[3,:] + ev[2,:])
                            ))

                gradPhi[3] = 0.
               # add the normal coupling term to get GradPhi
                #this fake gradPhi += vec(2*Cs*sum( sum(ev .* edge_normal_sum,dims=2) ./ sum(ev .* ev,dims=2) .* ev, dims=1))
               #gradPhi += 2*Cs*sum( sum(ev .* edge_normal_sum,dims=1) ./ sum(ev .* ev,dims=1) .* ev, dims=2)
               gradPhi[1] += 2*Cs*sum(filter(x -> !isnan(x), sum(ev[1,:] .* edge_normal_sum[1,:]) ./ sum(ev[1,:] .* ev[1,:]) .* ev[1,:] ))
               gradPhi[2] += 2*Cs*sum(filter(x -> !isnan(x), sum(ev[2,:] .* edge_normal_sum[2,:]) ./ sum(ev[2,:] .* ev[2,:]) .* ev[2,:] ))
               gradPhi[3] += 2*Cs*sum(filter(x -> !isnan(x), sum(ev[3,:] .* edge_normal_sum[3,:]) ./ sum(ev[3,:] .* ev[3,:]) .* ev[3,:] ))

               #gradPhi = gradPhi + ...
               #  2*Cs*sum( bsxfun(@times, sum(ev .* edge_normal_sum,2) ./ sum(ev .* ev,2) , ev), 1);

               # convert gradPhi to global coordinates
               #gradPhi = rotate(gradPhi, normals(i,:), 'to global');
               gradPhi = to_global(gradPhi,normals[:,i])
               # project to tangential plane
               # Matrix(1.0I, 3, 3) <- 3x3 Identity matrix
               tang_proj = Matrix(1.0I, 3, 3) - normals[:,i] * normals[:,i]'
               gradPhi = tang_proj * gradPhi

               P = normals[:,i] - 0.05*gradPhi
               normals[:,i] = P ./ norm(P)

               # println("edge number = ", i)
               # println("CDE = ", CDE[:,i])
               # println("gradPhi = ",gradPhi)
               #println("gradPhinorm = ",norm(gradPhi))
               #println("Phi = ",norm(gradPhi))
               #readline(stdin)
               #println(norm(gradPhi))
               #readline()
               if norm(gradPhi) < eps_inner
                   break
               end
               if k == max_iters_inner
                   println("inner iterations not converged")
               end
            end
        end
        #println("outer iters:")
        #println(maximum(sqrt.(sum(x -> x^2, normalsp - normals, dims=1))))

        print(" iter $m cost: ", maximum(sqrt.(sum(x -> x^2, normalsp - normals, dims=1))))
        if maximum(sqrt.(sum(x -> x^2, normalsp - normals, dims=1))) < eps_outer
            # biggest absolute change in normal vector
            println("paraboloid fit converged")
            break
        end

        if m==max_iters_outer
            println("WARNING!!! paraboloid fit not converged")
            println("WARNING!!! paraboloid fit not converged")
            println("WARNING!!! paraboloid fit not converged")
            println("WARNING!!! paraboloid fit not converged")
            println("WARNING!!! paraboloid fit not converged")
            println("WARNING!!! paraboloid fit not converged")
            println("WARNING!!! paraboloid fit not converged")
        end
    end
    println()
    return normals, CDE
end



# points, faces = expand_icosamesh(R=1, depth=2)
#
# #@load "./meshes/faces_critical_hyst_2_21.jld2" faces
# points = Array{Float64}(points)
# faces = Array{Int64}(faces)
# edges = make_edges(faces)
# connectivity = make_connectivity(edges)
# normals = Normals(points, faces)

(normalsd, CDEd) = make_deep_normals_spline(points, connectivity, edges, normals)
## testing deep normals spline
uugly_points = copy(points)
for i = 1:size(points,2)
    uugly_points[:,i]  /= norm(uugly_points[:,i])
end

n0 = Normals(uugly_points,faces)
@time n_old, CDE_old = make_normals_spline(uugly_points, connectivity, edges, n0)
@time n_new, CDE_new = make_deep_normals_spline(uugly_points, connectivity, edges, n0)
@time n_par, CDE_par, AB = make_normals_parab(uugly_points, connectivity, edges, n0; eps = 10^-8)

k1_old, k2_old = make_pc(CDE_old) # 1 theoretically
k1_new, k2_new = make_pc(CDE_new) # 1 theoretically
k1_par, k2_par = make_pc(CDE_par) # 1 theoretically

Plots.plot(k1_old.-1)
Plots.plot!(k1_new.-1)
Plots.plot!(k1_par.-1)

Plots.plot(k2_old.-1)
Plots.plot!(k2_new.-1)
Plots.plot!(k2_par.-1)


nerr_old = sqrt.(sum((n_old - uugly_points).^2,dims=1))
nerr_new = sqrt.(sum((n_new - uugly_points).^2,dims=1))
nerr_par = sqrt.(sum((n_par - uugly_points).^2,dims=1))

Plots.plot(nerr_old')
Plots.plot!(nerr_new')
Plots.plot!(nerr_par')


## uncoupled paraboloids
function make_normals_parab(points, connectivity, edges, normals0; eps = 10^-8)
    # described in Zinchenko et al. 1997
    #normals and CDE - parameters from locally fitted paraboloids
    #z = A*x+B*y+C*x^2+D*x*y+E*y^2
    #A and B should be 0 at the end
    println("fitting decoupled locals paraboloids")

    CDE = zeros(Float64, size(points)) # 3xN array
    AB = zeros(Float64, (2,size(points,2))) # 2xN array
    normals = copy(normals0)
    for i = 1:size(points,2)
        # get edges close to vertex
        # ev = 3xv matrix containing edge vectors from i-th to all adjecent points
        # ev is padded with 0-os
        #v = size(connectivity,1) ## maximum valence
        #ev = zeros(Float64,3,v)
        n0 = [0.,0.,0.]
        A, B, C, D, E = 0., 0., 0., 0., 0.
        while norm(n0 - normals[:,i]) > eps
            n0 = normals[:,i]
            b = zeros(Float64,(5,1))
            M = zeros(Float64,(5,5))
            for j = connectivity[:,i]
                if j == 0
                    break
                end

                ev = points[:,j] - points[:,i]
                ev = to_local(ev,normals[:,i])
                x,y,z = ev

                b += 2/(x^2+y^2+z^2) * (-z) * [x, y, x^2, x*y, y^2]
                M += 2/(x^2+y^2+z^2) * [x, y, x^2, x*y, y^2] * [x, y, x^2, x*y, y^2]'
            end
            A,B,C,D,E = M\(-b)
            normals[:,i] = to_global([ -A, -B, 1.] / sqrt(1. + A^2 + B^2), normals[:,i])

            # println(norm(n0 - normals[:,i]), "  ", i)
            # println(dot(normals[:,i] , points[:,i]))
            # println("A = ",A, " B = ", B)
            # println(C, " ", D, " ", E)
            # readline()
        end
        AB[:,i] = [A,B]
        CDE[:,i] = [C,D,E]
    end
    println("paraboloids fitted")
    return normals, CDE, AB
end

points, faces = expand_icosamesh(R=1, depth=2)

#@load "./meshes/faces_critical_hyst_2_21.jld2" faces
points = Array{Float64}(points)
faces = Array{Int64}(faces)

edges = make_edges(faces)
connectivity = make_connectivity(edges)
normals = Normals(points, faces)

@time n_par, CDE_par, AB = make_normals_parab(points, connectivity, edges, normals; eps = 10^-8)
