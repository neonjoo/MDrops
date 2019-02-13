function to_local(r::Array{Float64,1},normal::Array{Float64,1})
    cosf = normal[2] / sqrt( normal[1]^2 + normal[2]^2 )
    cost = normal[3]
    sinf = normal[1] / sqrt( normal[1]^2 + normal[2]^2 )
    sint = sqrt( 1 - normal[3]^2 )

    A = [cosf  -sinf  0;
        sinf*cost  cosf*cost  -sint;
        sinf*sint  cosf*sint  cost]

    rprim = A * r
    return rprim
end

function to_global(rprim::Array{Float64,1},normal::Array{Float64,1})
    cosf = normal[2] / sqrt( normal[1]^2 + normal[2]^2 )
    cost = normal[3]
    sinf = normal[1] / sqrt( normal[1]^2 + normal[2]^2 )
    sint = sqrt( 1 - normal[3]^2 )

    A = [cosf  -sinf  0;
        sinf*cost  cosf*cost  -sint;
        sinf*sint  cosf*sint  cost]

    r = A' * rprim
    return r
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

function make_connectivity(edges)
# find adjescent vertices
    valence = maximum(StatsBase.counts(edges)) # max number of adjescent vertices
    nvertices = maximum(edges)
    connectivity = zeros(Int64, valence, number_of_vertices) # create empty array padded with zeros
    for vert = 1:nvertices
        inds = findall(x -> vert in x, edges)
        for j = 1:size(inds,1)
            # write the other value that is not vert
            connectivity[j,vert] = edges[inds[j][1]%2+1,inds[j][2]]
        end
    end
    return connectivity
end

function make_edge_lens(points,edges)
    edge_lens = Array{Float64}(undef, size(edges,2))
    for i = 1:size(edges,2)
        edge_vec = points[:,edges[1,i]] - points[:,edges[2,i]]
        edge_lens[i] = norm(edge_vec)
    end
    return edge_lens
end

function make_normals_spline(points, connectivity, edges, normals0;
                    Cs=1.0, eps_inner=1e-5, eps_outer=1e-4,
                    max_iters_inner=1000, max_iters_outer=1000)
    #returns improved normals and CDE - parameters of locally fitted paraboloid
    #z = C*x^2+D*x*y+E*y^2
    #Cs is a coupling parameter. Zinchencko(2000) sets it to 1
    #eps inner and outer are convergence criteria.
    #eps outer >> eps_inner

    CDE = zeros(Float64, size(points)) # 3xN array
    gradPhi = [0.,0.,0.]
    normals = copy(normals0)
    #outer iterations
    for m = 1:max_iters_outer
        println(m)
        normalsp = copy(normals) #previous
        for i = 1:size(points,2)
            #inner iterations
            for k = 1:max_iters_inner
                # get edges close to vertex
                # ev = 3xv matrix containing edge vectors from i-th to all adjecent points
                # ev is padded with 0-os, where adjescent points < max adjescent points = v
                ev = zeros(Float64,3,size(connectivity,1))
                edge_normal_sum = zeros(Float64,3,size(connectivity,1))
                for (ind, j) in enumerate(connectivity[:,1])
                    # iterate through close points j
                    if j == 0
                        break
                    end
                    ev[:,ind] = points[:,j] - points[:,i]
                    edge_normal_sum[:,ind] = normals[:,j] + normals[:,i]

                    ev[:,ind] = to_local(ev[:,ind],normals[:,i])
                    edge_normal_sum[:,ind] = to_local(edge_normal_sum[:,ind],normals[:,i])
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

                # M * CDE  + b = 0;  M*CDE = -b
                #println(b)
                #readline(stdin)

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
               #println(norm(gradPhi))
               if norm(gradPhi) < eps_inner
                   break
               end
            end
        end
        println("outer iters:")
        println(maximum(sqrt.(sum(x -> x^2, normalsp - normals, dims=1))))
        if maximum(sqrt.(sum(x -> x^2, normalsp - normals, dims=1))) < eps_outer
            # biggest absolute change in normal vector
            break
        end
    end
    return normals, CDE
end
