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


function flip_edges!(faces, vertices, connectivity)
    maxx = size(vertices, 2)

    for i in 1:maxx
        i_num = length(filter(x -> x>0, connectivity[:,i]))
        if i_num <= 5
            #println("skipped $i-vertex: <=5 neighs")
            continue
        end

        for j in i+1:maxx
            if !(j in connectivity[:,i])
                #println("$i -/- $j !")
                continue
            end

            j_num = length(filter(x -> x>0, connectivity[:,j]))
            if j_num <= 5
                #println("vertex $i: skipped $j-vertex: <=5 neighs")
                continue
            end


            xi, xj = vertices[:,i], vertices[:,j]

            common = intersect(connectivity[:,i], connectivity[:,j])
            common = filter(x -> x>0, common)
            if length(common) != 2
                #println("skipping $i--$j, common len not 2, common = $common")
                continue
            end
            #println("$i--$j common = $common")
            k, m = common[1], common[2]
            xk, xm = vertices[:,k], vertices[:,m]

            kc = find_circumcenter(xi, xj, xk)
            mc = find_circumcenter(xi, xj, xm)

            # a threshold value
            d = norm(dot(xk-kc, xm-xk)) + norm(dot(xm-mc, xm-xk))

            if norm(xk - xm)^2 < d
                println("------------------------------------------------------------------- flipped $i--$j to $k--$m")

                readline(stdin)
                flip_connectivity!(faces, connectivity, i, j, k, m)
            end

        end # end j for

    end # end i for

end # end function


function flip_connectivity!(faces, connectivity, i, j, k, m)
    # adjusts faces & connectivity to the i--j  ->  k--m edge flip

    for s in 1:size(faces,2)
        # finds first(and only) column where all 3 indices appear and adjusts indices
        println("$s: $i $j $k // $(faces[:,s])")
        readline(stdin)

        if length(intersect([i,j,k], faces[:,s])) == 3
            row = findfirst(x-> x==j, faces[:,s])
            faces[row, s] = m

        elseif length(intersect([i,m,j], faces[:,s])) == 3
            row = findfirst(x-> x==i, faces[:,s])
            faces[row, s] = k
        end

    end # end s for

    for s in 1:size(connectivity,2)
        if length(intersect([j,k,m], connectivity[:,s])) == 3
            println("updatalino le conn table")
            row_j_in_i = findfirst(x-> x==j, connectivity[:,s])
            row_i_in_j = findfirst(x-> x==i, connectivity[:,j])
            row_k_in_m = findfirst(x-> x==0, connectivity[:,m])
            row_m_in_k = findfirst(x-> x==0, connectivity[:,k])

            connectivity[row_j_in_i, i] = 0
            connectivity[row_i_in_j, j] = 0

            if row_k_in_m == nothing || row_m_in_k == nothing
                println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                connectivity = vcat(connectivity, zeros(1, size(faces,2)))
                connectivity[end, m] = k
                connectivity[end, k] = m
            end

            connectivity[row_k_in_m, m] = k
            connectivity[row_m_in_k, k] = m

        end

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
