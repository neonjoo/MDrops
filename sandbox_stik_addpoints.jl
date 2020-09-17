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

datadir="/home/andris/sim_data/elongation_Bm5_lamdba10_mu30_adaptive_dt/"

files = readdir(datadir)
N = 25
file = files[3+N]
#println(file)
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
scatter(zface,critS)
scatter!(zface,critL)

function make_neighbor_faces(faces)
    neighbor_faces = Array{Int64}(undef, 3, 0)
    for i = 1:size(faces,2)

        vind_1, vind_2, vind_3 = faces[:,i] # vertex indices
        edge1, edge2, edge3 = sort([vind_1,vind_2]), sort([vind_1,vind_3]), sort([vind_2,vind_3]) #edges
        edges_i = [edge1,edge2,edge3]
        # find neighbors for triangle i
        this_neighbors = [0,0,0]
        this_neighbors_ind = 1
        for j = 1:size(faces,2)
            if i != j
                vjind_1, vjind_2, vjind_3 = faces[:,j] # vertex indices
                edgej1, edgej2, edgej3 = sort([vjind_1,vjind_2]), sort([vjind_1,vjind_3]), sort([vjind_2,vjind_3]) #edges

                # if any of the edges of j are in i
                if (edgej1 in edges_i) || (edgej2 in edges_i) || (edgej3 in edges_i)
                    this_neighbors[this_neighbors_ind] = j
                    this_neighbors_ind += 1
                    if this_neighbors_ind == 5
                        #this_neighbors_ind = 1
                        println("hmm, somtheing wrong - 4 neighbors")
                    end
                end
            end
        end

        neighbor_faces = cat(neighbor_faces,this_neighbors,dims=2)
    end
    return neighbor_faces
end

function mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = 0.55)
    k1, k2 = make_pc(CDE) # principal curvatures on vertices
    H = sqrt.(k1.^2 + k2.^2)

    ## mark faces that are too curved and need to be split
    marked_faces = falses(size(faces,2))

    for i = 1:size(faces,2)
        v1, v2, v3 = points[:,faces[1,i]], points[:,faces[2,i]], points[:,faces[3,i]] # triangle vertices
        d1, d2, d3 = norm(v1-v2), norm(v1-v3), norm(v2-v3) # edge lengths

        maxd = max(d1,d2,d3) # longest edge
        Hface = sum(H[faces[:,i]]) / 3 # average curvature H of vertices

        crit = Hface * maxd
        if crit > cutoff_crit
            marked_faces[i] = true
        end
    end

    # mark also all the faces that have at least two nearby marked faces
    neighbors_marked = false
    while !neighbors_marked
        marked_faces_old = copy(marked_faces)
        for i = 1:size(faces,2)
            if marked_faces[i] == false

                nearby_marked_faces = 0
                for j in neighbor_faces[:,i]
                    if marked_faces[j] == true
                        nearby_marked_faces += 1
                    end
                end
                if nearby_marked_faces > 1
                    marked_faces[i] == true
                end

            end
        end

        if marked_faces_old == marked_faces
            neighbors_marked = true
        end
    end # end while


    return marked_faces
end

neighbor_faces = make_neighbor_faces(faces)



function add_points(points, faces,normals, edges, CDE; cutoff_crit = 0.55)
    # add points in places where the curvature is too great
    println("started adding points")
    neighbor_faces = make_neighbor_faces(faces)
    marked_faces = mark_faces_for_splitting(points, faces, edges, CDE, neighbor_faces; cutoff_crit = cutoff_crit)
    println("Number of marked triangles: ",sum(marked_faces))

    faces_vec = [faces[:,i] for i in 1:size(faces,2)] # makes it easier to delete or add a face
    new_points = copy(points)
    checked_sides = Array{Int64}[] # create empty array
    is = Array{Int64}(undef,(1,3))
    js = Array{Int64}(undef,(1,3))
    ks = Array{Int64}(undef,(1,3))
    original_vertices_length = size(points,2)

    for i in 1:size(faces,2)
    if marked_faces[i] == true # loop throough face to be split
        #println("splitting marked triangle #$i")
        triangle = faces[:,i]

        is[1] = triangle[1]
        is[2] = triangle[2]
        is[3] = triangle[3]

        sides = [
            sort([triangle[1], triangle[2]]),
            sort([triangle[2], triangle[3]]),
            sort([triangle[3], triangle[1]])
        ]

        for j = 1:3
            if sides[j] in checked_sides
                js[j] = findfirst(x -> x==sides[j],checked_sides) + original_vertices_length
            else
                new_vertex = 0.5*( new_points[:,sides[j][1]] + new_points[:,sides[j][2]] )

                proj1 = project_on_given_paraboloid(CDE[:,sides[j][1]],normals[:,sides[j][1]], new_vertex, new_points[:,sides[j][1]])
                proj2 = project_on_given_paraboloid(CDE[:,sides[j][2]],normals[:,sides[j][2]], new_vertex, new_points[:,sides[j][2]])

                new_vertex = (proj1 + proj2)/2

                push!(checked_sides, sides[j]) # add side to checked sides
                new_points = hcat(new_points,new_vertex) # add to new vertices
                js[j] = size(checked_sides,1) + original_vertices_length
            end
        end

        # remove marked triangle from triangle array
        filter!(x->x!=triangle,faces_vec)

        # add triangles that the marked triangle was split into
        push!(faces_vec, [is[1],js[1],js[3]])
        push!(faces_vec, [js[1],is[2],js[2]])
        push!(faces_vec, [js[2],is[3],js[3]])
        push!(faces_vec, [js[1],js[2],js[3]])

        # loop through neighbor triangles and check if they are also marked
        # if not, split them in two

        for k in neighbor_faces[:,i]
        if marked_faces[k] == false
            # identify the point across the marked triangle
            neigh_side = [] # if neigh_point in triangle
            neigh_triangle = faces[:,k]
            for (neigh_ind, neigh_point) in enumerate(faces[:,k])
                if neigh_point in triangle
                    push!(neigh_side,neigh_point)
                else
                    #ks = neigh_point
                    # shift the indices in the neighboring triangle so that
                    # the point across the marked triangle is first
                    neigh_triangle = circshift(faces[:,k],4-neigh_ind)
                end
            end
            sort!(neigh_side)
            # point split in the middle of shared vertex
            jjs = findfirst(x -> x==neigh_side,checked_sides) + original_vertices_length
            # rename points in this triangle for ease of writing
            ks = neigh_triangle


            # remove adjacent triangle from triangle array
            filter!(x->x!=faces[:,k],faces_vec)
            # split it in two
            push!(faces_vec, [ks[2],jjs,ks[1]])
            push!(faces_vec, [jjs,ks[3],ks[1]])

        end
        end # k neighbor triangle loop

    end
    end # i triangle loop

    new_faces = []
    for (ind,face) in enumerate(faces_vec)
        if ind == 1
            new_faces = face
        else
            new_faces = hcat(new_faces,face)
        end
    end

    return new_points, new_faces
end

new_points, new_faces = add_points(points, faces,normals, edges, CDE; cutoff_crit = 0.55)

new_normals = Normals(new_points,new_faces)
