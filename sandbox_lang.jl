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


function test_flip(i::Int64, j::Int64, vertices, faces, connectivity)
    # checks whether edge between i,j vertices should be flipped
    # according to Zinchenko (2013)

    # number of neighbors
    i_num = length(filter(x -> x>0, connectivity[:,i]))
    j_num = length(filter(x -> x>0, connectivity[:,j]))

    # whether after flip each would have at least 5 neighbors
    if i_num <= 5 || j_num <= 5
        return false
    end

    xi, xj = vertices[:,i], vertices[:,j]

    common = intersect(connectivity[:,i], connectivity[:,j])
    common = filter(x -> x>0, common)

    k, m = common[1], common[2]
    xk, xm = vertices[:,k], vertices[:,m]

    kc = find_circumcenter(xi, xj, xk)
    mc = find_circumcenter(xi, xj, xm)

    # a threshold value
    d = norm(dot(xk-kc, xm-xk)) + norm(dot(xm-mc, xm-xk))

    if norm(xk - xm)^2 < d
        return true
    end

end


function flip_connectivity!(i, j, k, m, faces)
    # adjusts faces to the i-j  ->  k-m edge flip
        
    for s in 1:size(faces,2)
        if length(intersect([i,j,k], faces[:,s])) == 3
            row = findfirst(x-> x==j, faces[:,s])
            faces[row, s] = m

        elseif length(intersect([i,m,j], faces[:,s])) == 3
            row = findfirst(x-> x==i, faces[:,s])
            faces[row, s] = k
        end
    end

    return faces

end


function find_circumcenter(x1, x2, x3)
    # finds circumcenter coordinates of a triangle
    # according to wiki

    d = 2 * norm(cross(x1-x2, x2-x3))^2
    a = norm(x2-x3)^2 * dot(x1-x2, x1-x3) / d
    b = norm(x1-x3)^2 * dot(x2-x1, x2-x3) / d
    c = norm(x1-x2)^2 * dot(x3-x1, x3-x2) / d

    return a*x1 + b*x2 + c*x3

end
