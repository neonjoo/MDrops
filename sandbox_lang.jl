function make_connectivity(edges)
# find adjescent vertices
    valence = maximum(StatsBase.counts(edges)) # max number of adjescent vertices
    nvertices = maximum(edges)
    connectivity = zeros(Int64, valence, nvertices) # create empty array padded with zeros
    for vert = 1:nvertices
        inds = findall(x -> vert in x, edges)
        for j = 1:size(inds,1)
            # write the other value that is not vert
            connectivity[j,vert] = edges[inds[j][1]%2+1,inds[j][2]]
        end
    end
    return connectivity
end



function make_edges(faces)
    edges = Array{Int64}(undef, 2, 0)
    for i in 1:size(faces,2) # faces
        face = faces[:,i]
        for j in 1:size(face,1) # vertices in face
            edge = sort([face[j],face[j%3+1]])

            ##check if edge is in edges
            duplicate = false
            for k = 1:size(edges,2)
                if edge == edges[:,k]
                    duplicate = true
                    break
                end
            end

            if !duplicate
                edges = cat(edges,edge,dims=2)
            end

        end
    end
    return edges
end


function flip_edges!(faces, connectivity, vertices)
    # flip edges to improve mesh, updates "faces" and "connectivity"

    maxx = size(vertices, 2)
    global continue_flip = true

    while continue_flip
        global continue_flip
        continue_flip = false
        println("started flip loop")

        for i in 1:maxx
            # num of i-th vertex neighbors
            i_num = length(filter(x -> x>0, connectivity[:,i]))
            if i_num <= 5
                continue
            end

            for j in i+1:maxx
                if !(j in connectivity[:,i])
                    continue
                end

                # num of j-th vertex neighbors
                j_num = length(filter(x -> x>0, connectivity[:,j]))
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
                    println("--------------------- flippening $i--$j to $k--$m")
                    #readline(stdin)
                    flip_connectivity!(faces, connectivity, i, j, k, m)
                    continue_flip = true
                    # break
                end

            end # end j for
            if continue_flip
                break
            end
        end # end i for

    end # end while

end # end function


function flip_connectivity!(faces, connectivity, i, j, k, m)
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


    for s in 1:size(connectivity,2)
        if length(intersect([j,k,m], connectivity[:,s])) == 3
            # row of j in i-th column etc.
            row_j_in_i = findfirst(x-> x==j, connectivity[:,s])
            row_i_in_j = findfirst(x-> x==i, connectivity[:,j])
            row_k_in_m = findfirst(x-> x==0, connectivity[:,m])
            row_m_in_k = findfirst(x-> x==0, connectivity[:,k])

            # cut the i--j edge in "connectivity" by sliding column values up 1 row, placing 0 in the end
            connectivity[row_i_in_j:end, j] = [connectivity[row_i_in_j+1:end, j]; 0]
            connectivity[row_j_in_i:end, s] = [connectivity[row_j_in_i+1:end, s]; 0]

            # adds a zero row if either of new vertices k or m are connected to some previous maximum connectivity (aka valence)
            if row_k_in_m == nothing || row_m_in_k == nothing
                println("padded row of 0s on bottom of \"connectivity\"")
                connectivity = vcat(connectivity, zeros(1, size(connectivity,2)))
                connectivity[end, m] = k
                connectivity[end, k] = m
                continue
            end

            connectivity[row_k_in_m, m] = k
            connectivity[row_m_in_k, k] = m

        end # end if

    end # end s for




end # end function


function find_circumcenter(x1, x2, x3)
    # finds circumcenter coordinates of a triangle
    # according to wiki

    d = 2 * norm(cross(x1-x2, x2-x3))^2
    a = norm(x2-x3)^2 * dot(x1-x2, x1-x3) / d
    b = norm(x1-x3)^2 * dot(x2-x1, x2-x3) / d
    c = norm(x1-x2)^2 * dot(x3-x1, x3-x2) / d

    return a*x1 + b*x2 + c*x3

end
