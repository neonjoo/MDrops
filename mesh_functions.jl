function Normals(points, faces)
    # return normals at vertices as weighted average of surrounding triangle normals
    # triangle angles are the weights
    all_normals = zeros(3, size(points,2))
    for i in 1:size(points)[2]
        mask = findall(x-> x == i, faces)
        neighbors = zeros(Int64, 3, length(mask))
        num_neighbors = length(mask)
        for n in 1:num_neighbors
            neighbors[:, n] = faces[:,mask[n][2]]
        end
        normal_i = zeros(3, 1)
        for j in 1:num_neighbors
            all_i = findall(neighbors .== i)
            v1 = neighbors[all_i[j][1], all_i[j][2]]
            v2 = neighbors[(all_i[j][1]) % 3 + 1, all_i[j][2]]
            v3 = neighbors[(all_i[j][1] + 1) % 3 + 1, all_i[j][2]]
            #println("v1 = $(norm(v1)), v2 = $(norm(v2)), v3 = $(norm(v3))")
            vec1 = points[:, v3] - points[:, v1]
            vec2 = points[:, v2] - points[:, v1]
            #println(norm(vec1), norm(vec2))
            normal_j = cross(vec1, vec2)
            angle_j = acos(dot(vec1, vec2) / norm(vec1) / norm(vec2))
            normal_i += angle_j * normal_j
        end
        normal_i = normalize(normal_i[:,1])
        all_normals[:, i] = normal_i
    end
    return all_normals
end

function make_icosamesh(;R=1)
    vertices = Array{Float64}(undef, (12,3))
    triangles = Array{Int64}(undef, (20,3))

    #R = 1
    # x,y,z coordinates of the vertices
    vertices[1,:] = [0,0,R]
    vertices[12,:] = [0,0,-R]

    ϕ = 0.
    for i in 2:11

        if i%2 == 0
            θ = 2*asin(1/(2*sin(2pi/5)))
        elseif i%2 == 1
            θ = pi - 2*asin(1/(2*sin(2pi/5)))
        end

        vertices[i,:] = R*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ) ]

        ϕ += pi/5
        #global ϕ += pi/5
    end

    #triangles are marked counter clockwise, looking from the outside
    #top
    triangles[1,:] = [1,2,10]
    triangles[2,:] = [1,4,2]
    triangles[3,:] = [1,6,4]
    triangles[4,:] = [1,8,6]
    triangles[5,:] = [1,10,8]
    #middle
    triangles[6,:] = [2,3,11]
    triangles[7,:] = [2,4,3]
    triangles[8,:] = [3,4,5]
    triangles[9,:] = [4,6,5]
    triangles[10,:] = [5,6,7]
    triangles[11,:] = [6,8,7]
    triangles[12,:] = [7,8,9]
    triangles[13,:] = [8,10,9]
    triangles[14,:] = [9,10,11]
    triangles[15,:] = [2,11,10]
    #bottom
    triangles[16,:] = [3,12,11]
    triangles[17,:] = [3,5,12]
    triangles[18,:] = [5,7,12]
    triangles[19,:] = [7,9,12]
    triangles[20,:] = [9,11,12]

    return vertices, triangles
end

function expand_icosamesh(;R=1,depth=3)
    vertices, triangles = make_icosamesh(R=R)
    for iter = 1:depth
        checked_sides = Array{Int64}[] # create empty array
        new_triangles = Array{Int64}(undef,(0,3)) # create 0x3 uninitialized matrix
        is = Array{Int64}(undef,(1,3))
        js = Array{Int64}(undef,(1,3))
        original_vertices_length = size(vertices,1)
        for i = 1:size(triangles,1)
            triangle = triangles[i,:]

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
                    new_vertex = 0.5*( vertices[sides[j][1],:] + vertices[sides[j][2],:] )
                    scale = R/sqrt(new_vertex[1]^2 + new_vertex[2]^2 + new_vertex[3]^2)
                    new_vertex *= scale
                    push!(checked_sides, sides[j]) # add side to checked sides
                    vertices = vcat(vertices,new_vertex') # add to new vertices
                    js[j] = size(checked_sides,1) + original_vertices_length
                end
            end

            new_triangles = vcat(new_triangles,[is[1],js[1],js[3]]')
            new_triangles = vcat(new_triangles,[js[1],is[2],js[2]]')
            new_triangles = vcat(new_triangles,[js[2],is[3],js[3]]')
            new_triangles = vcat(new_triangles,[js[1],js[2],js[3]]')
        end
        triangles = new_triangles
    end
    return vertices', triangles'
end

function make_pc(CDE::Array{Float64,2})
    # returns principal curvatures k1, k2 at the vertices
    # with locally fitted paraboloids z = Cx^2 + Dxy + Ey^2
    # k1 > k2

    C = CDE[1,:]
    D = CDE[2,:]
    E = CDE[3,:]

    k1 = -C - E + sqrt.((C-E).^2 + D.^2)
    k2 = -C - E - sqrt.((C-E).^2 + D.^2)

    # the below formulas gave errors small negative numbers under sqrt
    #k1 = -C - E + sqrt.(C.^2 + D.^2 - 2*C.*E + E.^2)
    #k2 = -C - E - sqrt.(C.^2 + D.^2 - 2*C.*E + E.^2)

    return k1,k2
end

function make_dS(points,faces)
    dS = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        dS[v1] += area/3
        dS[v2] += area/3
        dS[v3] += area/3
    end
    return dS
end

function make_S(points,faces)
    S = 0
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        S += area
    end
    return S
end

function make_pc_local(CDE_local::Array{Float64,1},x::Float64,y::Float64)
    # returns principal curvatures k1, k2 at the point ( x , y , z(x,y) )
    # on a locally fitted paraboloid z = Cx^2 + Dxy + Ey^2
    # k1 > k2 (hopefully)

    C = CDE_local[1]
    D = CDE_local[2]
    E = CDE_local[3]

    magN = 1 + (2*C*x + D*y)^2 + (D*x + 2*E*y)^2; # a repeating value

    k2 = -1/magN^2 *
        (
        C*sqrt(magN) + E*sqrt(magN) - C*D^2*x^2*sqrt(magN) + 4*C^2*E*x^2*sqrt(magN) -
        D^3*x*y*sqrt(magN) + 4*C*D*E*x*y*sqrt(magN) - D^2*E*y^2*sqrt(magN) + 4*C*E^2*y^2*sqrt(magN) +
        0.5*sqrt(abs(
            4*(D^2 - 4*C*E)*(1 + 4*C^2*x^2 + 4*C*D*x*y + 4*D*E*x*y + 4*E^2*y^2 + D^2*(x^2 + y^2))^2 +
            4*(1 + (2*C*x + D*y)^2 + (D*x + 2*E*y)^2) * (E + 4*C^2*E*x^2 - D^3*x*y - D^2*E*y^2 + C*(1 - D^2*x^2 + 4*D*E*x*y + 4*E^2*y^2))^2
            ))
        )

    k1 = -1/magN^2 *
        (
        C*sqrt(magN) + E*sqrt(magN) - C*D^2*x^2*sqrt(magN) + 4*C^2*E*x^2*sqrt(magN) -
        D^3*x*y*sqrt(magN) + 4*C*D*E*x*y*sqrt(magN) - D^2*E*y^2*sqrt(magN) + 4*C*E^2*y^2*sqrt(magN) -
        0.5*sqrt(abs(
            4*(D^2 - 4*C*E)*(1 + 4*C^2*x^2 + 4*C*D*x*y + 4*D*E*x*y + 4*E^2*y^2 + D^2*(x^2 + y^2))^2 +
            4*(1 + (2*C*x + D*y)^2 + (D*x + 2*E*y)^2) * (E + 4*C^2*E*x^2 - D^3*x*y - D^2*E*y^2 + C*(1 - D^2*x^2 + 4*D*E*x*y + 4*E^2*y^2))^2
            ))
        )

    return k1,k2
end

function to_local(r::Array{Float64,1},normal::Array{Float64,1})
    # rotate a vector to local coordinate system
    # with z axis along a normal
    cosf = normal[2] / sqrt( normal[1]^2 + normal[2]^2)
    cost = normal[3]
    sinf = normal[1] / sqrt( normal[1]^2 + normal[2]^2)
    sint = sqrt( 1 - normal[3]^2 )

    A = [cosf  -sinf  0;
        sinf*cost  cosf*cost  -sint;
        sinf*sint  cosf*sint  cost]

    rprim = A * r
    return rprim
end

function to_global(rprim::Array{Float64,1},normal::Array{Float64,1})
    # rotate a vector back to global coordinates
    # from a local coordinate system along a normal
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

function make_closefaces(faces)
    # find adjescent triangles
    valence = maximum(StatsBase.counts(faces)) # max number of adjescent vertices
    nvertices = maximum(faces)
    closefaces = zeros(Int64, valence, nvertices) # create empty array padded with zeros
    for vert = 1:nvertices
        inds = findall(x -> vert in x, faces)
        for j = 1:size(inds,1)
            # write the other value that is not vert
            closefaces[j,vert] = inds[j][2]
        end
    end
    return closefaces
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
                    max_iters_inner=1000, max_iters_outer=100)
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

                # get edges close to vertex
                # ev = 3xv matrix containing edge vectors from i-th to all adjecent points
                # ev is padded with 0-os, where adjescent points < max adjescent points = v
                ev = zeros(Float64,3,size(connectivity,1))
                edge_normal_sum = zeros(Float64,3,size(connectivity,1))
                for (ind, j) in enumerate(connectivity[:,i])
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

               if norm(gradPhi) < eps_inner
                   break
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

function make_min_edges(points,connectivity)
    # finds the length of the smallest of adjecscent edges to a point
    min_edges = ones(Float64,size(points,2)) * 1e200
    for i = 1:size(points,2)
        for j in connectivity[:,i]
            if j == 0
                break
            end
            #edge vector length
            evl = norm(points[:,j] - points[:,i])
            if min_edges[i] > evl
                min_edges[i] = evl
            end
        end
    end
    return min_edges
end

function project_on_drop(points::Array{Float64,2},CDE::Array{Float64,2},normals::Array{Float64,2},r0::Array{Float64,1})
    # projects a point r0 on the closet fitted paraboloid
    i = argmin(sum((points .- r0).^2,dims=1))[2]

    r = points[:,i]

    r0 = to_local(r0 - r, normals[:,i])

    f(x::Array{Float64,1}) = (x[1]-r0[1])^2 + (x[2]-r0[2])^2 +
            (CDE[1,i]*x[1]^2 + CDE[2,i]*x[1]*x[2] + CDE[3,i]*x[2]^2 - r0[3])^2
    function g!(storage::Array{Float64,1},x::Array{Float64,1}) # gradient of f
        grad = storage
        grad[1] = 2*(x[1]-r0[1] + (2*CDE[1,i]*x[1] + CDE[2,i]*x[2])*(CDE[1,i]*x[1]^2 + CDE[2,i]*x[1]*x[2] + CDE[3,i]*x[2]^2 - r0[3]))
        grad[2] = 2*(x[2]-r0[2] + (CDE[2,i]*x[1] + 2*CDE[3,i]*x[2])*(CDE[1,i]*x[1]^2 + CDE[2,i]*x[1]*x[2] + CDE[3,i]*x[2]^2 - r0[3]))
    end

    function h!(storage::Array{Float64,2},x::Array{Float64,1}) # hessian of f
        h = storage
        h[1, 1] = 2.0 + 2*(2*CDE[1,i]*x[1] + CDE[2,i]*x[2])^2 + 4*CDE[1,i]*(CDE[1,i]*x[1]^2 + CDE[2,i]*x[1]*x[2] + CDE[3,i]*x[2]^2 - r0[3])
        h[1, 2] = 2*CDE[1,i]*x[1]*(3*CDE[2,i]*x[1] + 4*CDE[3,i]*x[2]) + 2*CDE[2,i]*(2*CDE[2,i]*x[1]*x[2] + 3*CDE[3,i]*x[2]^2 - r0[3])
        h[2, 1] = 2*CDE[1,i]*x[1]*(3*CDE[2,i]*x[1] + 4*CDE[3,i]*x[2]) + 2*CDE[2,i]*(2*CDE[2,i]*x[1]*x[2] + 3*CDE[3,i]*x[2]^2 - r0[3])
        h[2, 2] = 2.0 + 2*(CDE[2,i]*x[1] + 2*CDE[3,i]*x[2])^2 + 4*CDE[3,i]*(CDE[1,i]*x[1]^2 + CDE[2,i]*x[1]*x[2] + CDE[3,i]*x[2]^2 - r0[3])
    end

    x0 = [r0[1],r0[2]]

    res = Optim.optimize(f,g!,h!,x0,NewtonTrustRegion(),Optim.Options(f_tol=1.e-8))
    #println(res)
    x1 = Optim.minimizer(res)[1]
    x2 = Optim.minimizer(res)[2]
    r0 = [x1,x2,CDE[1,i]*x1^2 + CDE[2,i]*x1*x2 + CDE[3,i]*x2^2]
    #println(r0)
    #println(x0)

    r0 = to_global(r0,normals[:,i]) + r
    return r0
end

function project_on_given_paraboloid(CDE::Array{Float64,1},normal::Array{Float64,1},r0::Array{Float64,1}, r)
    # projects a point r0 on the paraboloid fitted around point at r
    # z = Cx^2 + D*x*y + E*y^2
    r0 = to_local(r0 - r, normal)

    f(x::Array{Float64,1}) = (x[1]-r0[1])^2 + (x[2]-r0[2])^2 +
            (CDE[1]*x[1]^2 + CDE[2]*x[1]*x[2] + CDE[3]*x[2]^2 - r0[3])^2

    function g!(storage::Array{Float64,1},x::Array{Float64,1}) # gradient of f
        grad = storage
        grad[1] = 2*(x[1]-r0[1] + (2*CDE[1]*x[1] + CDE[2]*x[2])*(CDE[1]*x[1]^2 + CDE[2]*x[1]*x[2] + CDE[3]*x[2]^2 - r0[3]))
        grad[2] = 2*(x[2]-r0[2] + (CDE[2]*x[1] + 2*CDE[3]*x[2])*(CDE[1]*x[1]^2 + CDE[2]*x[1]*x[2] + CDE[3]*x[2]^2 - r0[3]))
    end

    function h!(storage::Array{Float64,2},x::Array{Float64,1}) # hessian of f
        h = storage
        h[1, 1] = 2.0 + 2*(2*CDE[1]*x[1] + CDE[2]*x[2])^2 + 4*CDE[1]*(CDE[1]*x[1]^2 + CDE[2]*x[1]*x[2] + CDE[3]*x[2]^2 - r0[3])
        h[1, 2] = 2*CDE[1]*x[1]*(3*CDE[2]*x[1] + 4*CDE[3]*x[2]) + 2*CDE[2]*(2*CDE[2]*x[1]*x[2] + 3*CDE[3]*x[2]^2 - r0[3])
        h[2, 1] = 2*CDE[1]*x[1]*(3*CDE[2]*x[1] + 4*CDE[3]*x[2]) + 2*CDE[2]*(2*CDE[2]*x[1]*x[2] + 3*CDE[3]*x[2]^2 - r0[3])
        h[2, 2] = 2.0 + 2*(CDE[2]*x[1] + 2*CDE[3]*x[2])^2 + 4*CDE[3]*(CDE[1]*x[1]^2 + CDE[2]*x[1]*x[2] + CDE[3]*x[2]^2 - r0[3])
    end

    x0 = [r0[1],r0[2]]

    res = Optim.optimize(f,g!,h!,x0,NewtonTrustRegion(),Optim.Options(f_tol=1.e-8))
    x1 = Optim.minimizer(res)[1]
    x2 = Optim.minimizer(res)[2]
    r0 = [x1,x2,CDE[1]*x1^2 + CDE[2]*x1*x2 + CDE[3]*x2^2]

    r0 = to_global(r0,normal) + r
    return r0
end

function active_stabilize(points0::Array{Float64,2},faces::Array{Int64,2},CDE::Array{Float64,2},connectivity::Array{Int64,2},edges::Array{Int64,2},normals::Array{Float64,2};
    deltakoef=0.01, R0=1.0, gamma=0.25, p=50, r=100, checkiters=100, maxiters=1000,critSc = 0.75,critCdelta = 1.15)
    # actively rearange vertices on a surfaces given by fitted paraboloids
    # as per Zinchenko(2013)
    println("active stabilization")
    points = copy(points0)

    closefaces = make_closefaces(faces)
    dS = make_dS(points,faces)
    k1,k2 = make_pc(CDE)
    LAMBDA = k1.^2 + k2.^2 .+ 0.004/R0^2
    K = 4/(sqrt(3) * size(faces,2)) * sum(LAMBDA.^gamma .* dS)
    hsq = K * LAMBDA.^(-gamma)

    no_improvement = true
    # initiallize improvement criteria
    xij = make_edge_lens(points,edges)
    hij = sqrt.(0.5*(hsq[edges[1,:]].^2 + hsq[edges[2,:]].^2))

    Sc0 = maximum(xij./hij) / minimum(xij./hij)
    Cdelta_min0 = minimum(make_Cdeltas(points, faces))
    for iter = 1:maxiters
        #println(iter)
        gradE = make_gradE(points,faces,closefaces,hsq; p=p,r=r)
        delta = deltakoef * minimum(make_min_edges(points,connectivity) ./ sum(sqrt.(gradE.^2),dims=1))

        #println("dPoints= ",delta*maximum(sum(sqrt.(gradE.^2),dims=1)))
        #println("E = ", make_E(points,faces,hsq; p=p,r=r))


        points = points - delta * gradE

        #project points on the drop
        for i = 1:size(points,2)
            points[:,i] = project_on_drop(points0,CDE,normals,points[:,i])
        end

        # recalculate the mesh parameters
        dS = make_dS(points,faces)
        #k1,k2 = make_pc(CDE)
        for i = 1:size(points,2)
            r0 = points[:,i]
            minind = argmin(sum((points0 .- r0).^2,dims=1))[2]
            r0 = to_local(r0-points0[:,minind],normals[:,minind])
            k1[i],k2[i] =  make_pc_local(CDE[:,minind],r0[1],r0[2])
        end
        LAMBDA = k1.^2 + k2.^2 .+ 0.004/R0^2
        K = 4/(sqrt(3) * size(faces,2)) * sum(LAMBDA.^gamma .* dS)
        hsq = K * LAMBDA.^(-gamma)
        if iter < checkiters
            xij = make_edge_lens(points,edges)
            hij = sqrt.(0.5*(hsq[edges[1,:]].^2 + hsq[edges[2,:]].^2))

            Sc = maximum(xij./hij) / minimum(xij./hij)
            Cdelta_min = minimum(make_Cdeltas(points, faces))

            #println("Sc/Sc0 = ",Sc/Sc0)
            #println("Cdelta/Cdelta0 = ",Cdelta_min/Cdelta_min0)
            if Sc > critSc*Sc0 && Cdelta_min < critCdelta*Cdelta_min0
                no_improvement = false
            end
        end

        if iter == checkiters
            if no_improvement == true
                println("no significant improvement achieved")
                println("reversing changes")
                points = points0
                break
            else
                println("improvement detected in the first ", checkiters, " iterations")
                println("iterating for ", maxiters - checkiters, " more iterations")
            end
        end

        # e = make_E(points, faces, hsq,p=p, r=r)
        # println("E = $e")

        if iter%500 == 0
            println("iteration ",iter)

        end
    end
    return points
end

function make_gradE(points,faces,closefaces,hsq; p=50,r=100)
    Cet = sqrt(3) / 12 # target Cdelta value
    Nvertices = size(points, 2)
    gradE = zeros(3,Nvertices)
    for i = 1:Nvertices
        for faceind in closefaces[:,i]
            if faceind == 0
                break
            end
        # order radius vectors x so that first one points to vertex i
        iind = findfirst(x->x==i,faces[:,faceind])
        x1 = points[:,faces[iind,faceind]]
        x2 = points[:,faces[iind%3+1,faceind]]
        x3 = points[:,faces[(iind+1)%3+1,faceind]]

        hsq1 = hsq[faces[iind,faceind]]
        hsq2 = hsq[faces[iind%3+1,faceind]]
        hsq3 = hsq[faces[(iind+1)%3+1,faceind]]

        a = norm(x2 - x1)
        b = norm(x3 - x2)
        c = norm(x1 - x3)

        Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        x12 = x2 - x1
        x13 = x3 - x1
        h12sq = 0.5*(hsq1 + hsq2)
        h13sq = 0.5*(hsq1 + hsq3)

        dCdx = 1/(4*Cdelta*(a^2 + b^2 + c^2)^3) * (
                -(a^4 + b^4 + c^4 - (a^2 + b^2 + c^2)*a^2) * x12
                -(a^4 + b^4 + c^4 - (a^2 + b^2 + c^2)*c^2) * x13
                )

        gradE[:,i] = gradE[:,i] +
            r*(Cet/Cdelta)^(r-1) * Cet * (-1/Cdelta^2) * dCdx +
            0.5*p*( 0.5*( dot(x12,x12)/h12sq + h12sq/dot(x12,x12)) )^(p-1) *(-1)* (1/h12sq - h12sq/dot(x12,x12)^2) * x12 +
            0.5*p*( 0.5*( dot(x13,x13)/h13sq + h13sq/dot(x13,x13)) )^(p-1) *(-1)* (1/h13sq - h13sq/dot(x13,x13)^2) * x13
        end
    end
    return gradE
end

function make_E(points,faces,hsq; p=50,r=100)
    Cet = sqrt(3) / 12 # target Cdelta value
    Ntriangles = size(faces, 2)

    E = 0
    for i = 1:Ntriangles
        x1 = points[:,faces[1,i]]
        x2 = points[:,faces[2,i]]
        x3 = points[:,faces[3,i]]

        hsq1 = hsq[faces[1,i]]
        hsq2 = hsq[faces[2,i]]
        hsq3 = hsq[faces[3,i]]

        a = norm(x2 - x1)
        b = norm(x3 - x2)
        c = norm(x1 - x3)

        hsqa = 0.5*( hsq2 + hsq1 )
        hsqb = 0.5*( hsq3 + hsq2 )
        hsqc = 0.5*( hsq1 + hsq3 )

        Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

        E += (Cet/Cdelta)^r +
            0.5 * (0.5 * (a^2/hsqa + hsqa/a^2) )^p +
            0.5 * (0.5 * (b^2/hsqb + hsqb/b^2) )^p +
            0.5 * (0.5 * (c^2/hsqc + hsqc/c^2) )^p

    end
    return E
end

function make_Cdeltas(points, faces)
    Nfaces = size(faces,2)
    Cdeltas = zeros(Nfaces)

    for i in 1:Nfaces


        a = norm(points[:,faces[2,i]] - points[:,faces[1,i]])
        b = norm(points[:,faces[3,i]] - points[:,faces[2,i]])
        c = norm(points[:,faces[1,i]] - points[:,faces[3,i]])

        Cdeltas[i] = 0.25*sqrt( 1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2 )
    end
    return Cdeltas
end

function flip_edges(faces, connectivity, vertices; max_neighbors=7)
    # flip edges to improve mesh, updates "faces" and "connectivity"

    println("flipping edges")
    maxx = size(vertices, 2)
    global continue_flip = true
    global flipped_any = false
    max_flips = 10
    flips_tried = 0
    while continue_flip & (flips_tried < max_flips)
        global continue_flip
        continue_flip = false
        flips_tried += 1


        for i in 1:maxx
            # num of i-th vertex neighbors
            i_num = length(filter(x -> x>0, connectivity[:,i]))
            if i_num <= 5
                continue
            end

            for j in filter(x -> x>0, connectivity[:,i])

                # num of j-th vertex neighbors
                j_num = length(filter(x -> x>0, connectivity[:,j]))
                if j_num <= 5
                    continue
                end

                xi, xj = vertices[:,i], vertices[:,j]

                common = intersect(connectivity[:,i], connectivity[:,j])
                common = filter(x -> x>0, common)

                if length(common) != 2
                    println("This should never print if the mesh is OK")
                    continue
                end

                k, m = common[1], common[2]

                k_num = length(filter(x -> x>0, connectivity[:,k]))
                if k_num >= max_neighbors
                    continue
                end

                m_num = length(filter(x -> x>0, connectivity[:,m]))
                if m_num >= max_neighbors
                    continue
                end

                xk, xm = vertices[:,k], vertices[:,m]

                kc = find_circumcenter(xi, xj, xk)
                mc = find_circumcenter(xi, xj, xm)

                # a threshold value for flipping (Zinchenko2013)
                d = norm(dot(xk-kc, xm-xk)) + norm(dot(xm-mc, xm-xk))

                if norm(xk - xm)^2 < d
                    #println("--------------------- flippening $i--$j to $k--$m")
                    faces, connectivity = flip_connectivity(faces, connectivity, i, j, k, m)
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

    return faces, connectivity, flipped_any
end

function flip_connectivity(faces, connectivity, i, j, k, m)
    # adjusts faces & connectivity to the i--j  ->  k--m edge flip

    #println()
    #println("entered flip_con ---------------------")
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

    return faces, connectivity
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

function make_F(triangles, vertices, V)

    # triangles = triangles'
    # vertices = vertices'
    # V = V'

    Ntriangles = size(triangles, 2)

    F = 0

    for i = 1:Ntriangles
        x = [vertices[:,triangles[1,i]],
             vertices[:,triangles[2,i]],
             vertices[:,triangles[3,i]]]

        v = [V[:,triangles[1,i]],
             V[:,triangles[2,i]],
             V[:,triangles[3,i]]]

        # this_hsq = [hsq[triangles[i,1]]
        #             hsq[triangles[i,2]]
        #             hsq[triangles[i,3]]]

        a = norm(x[2] - x[1])
        b = norm(x[3] - x[2])
        c = norm(x[1] - x[3])
        # println("$i: $a $b $c")
        # try
        #     global Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        # catch
        #     println("sadness")
        #     global Cdelta = 0
        # end

        Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

        A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
            (4*Cdelta * (a^2 + b^2 + c^2)^3)
        B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
            (4*Cdelta * (a^2 + b^2 + c^2)^3)
        C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
            (4*Cdelta * (a^2 + b^2 + c^2)^3);

        dCdeltadt = -A * dot(x[2] - x[1], v[2] - v[1]) +
                    -B * dot(x[3] - x[2], v[3] - v[2]) +
                    -C * dot(x[1] - x[3], v[1] - v[3])


        # F = F + ...
        #     0.4 / Cdelta^2 * dCdeltadt^2 + ...
        #     2*( dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[2])) - 0.5*(this_hsq[1] + this_hsq[2]) / a^4) )^2 + ...
        #     2*( dot(x[3,:] - x[2,:], v[3,:] - v[2,:]) * (1/(0.5*(this_hsq[3] + this_hsq[2])) - 0.5*(this_hsq[3] + this_hsq[2]) / b^4) )^2 + ...
        #     2*( dot(x[1,:] - x[3,:], v[1,:] - v[3,:]) * (1/(0.5*(this_hsq[1] + this_hsq[3])) - 0.5*(this_hsq[1] + this_hsq[3]) / c^4) )^2

        F = F + 0.4 / Cdelta^2 * dCdeltadt^2 +
            2*( dot(x[2] - x[1], v[2] - v[1]))^2 +
            2*( dot(x[3] - x[2], v[3] - v[2]))^2 +
            2*( dot(x[1] - x[3], v[1] - v[3]))^2

    end

    return F
end

function make_tanggradF(normals,triangles, vertices, V)

    # normals = normals'
    # triangles = triangles'
    # vertices = vertices'
    # V = V'

    Nvertices = size(vertices, 2)
    gradF = zeros(3,Nvertices)

    for i = 1:Nvertices

        # finds the i-th triangle indices in the triangle 2D array
        mask = findall(x-> x == i, faces)
        num = length(mask)
        (row, col) = zeros(Int64, 1, num), zeros(Int64, 1, num)
        for n in 1:num
            row[n], col[n] = mask[n][1], mask[n][2]
        end


        for this_triang = 1:length(row)
            j1 = (row[this_triang])%(3) + 1 # from 1 to 3
            j2 = (row[this_triang] + 1)%(3) +1 # from 1 to 3


            # radiusvektori uz sho trijsturu virsotneem
            # pirmais atbilst i

            # println()
            # println("triang $i, kaims = $this_triang ----------------------------------------------------------------------------")
            #
            # println()
            # println("x index")
            # println("trig size = ", size(triangles))
            # println("vert size = ", size(vertices))
            #
            # println("$(row[this_triang]), $(col[this_triang]) -> ", triangles[row[this_triang], col[this_triang]])
            # println("$j1, $(col[this_triang]) -> ", triangles[j1, col[this_triang]])
            # println("$j2, $(col[this_triang]) -> ", triangles[j2, col[this_triang]])
            #
            # println(size(triangles))
            # #println(triangles'[1:20, :])
            # println(triangles[1, 1:15])
            # println(triangles[2, 1:15])
            # println(triangles[3, 1:15])
            # println()


            x = [vertices[:,triangles[row[this_triang],col[this_triang]]],
                vertices[:,triangles[j1, col[this_triang]]],
                vertices[:,triangles[j2, col[this_triang]]]]

            v = [V[:, triangles[row[this_triang],col[this_triang]]],
                V[:, triangles[j1, col[this_triang]]],
                V[:, triangles[j2, col[this_triang]]]]

            # this_hsq = [hsq[triangles[row[this_triang],col[this_triang]]]
            #     hsq[triangles[row[this_triang],j1]]
            #     hsq[triangles[row[this_triang],j2]]]

            a = norm(x[2] - x[1])
            b = norm(x[3] - x[2])
            c = norm(x[1] - x[3])

            Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

            A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)
            B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)
            C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)

            dCdeltadt = -A * dot(x[2] - x[1], v[2] - v[1]) +
                        -B * dot(x[3] - x[2], v[3] - v[2]) +
                        -C * dot(x[1] - x[3], v[1] - v[3])

            # gradF[i,:] = gradF[i,:] + ...
            #             0.4 / Cdelta^2 * 2 * dCdeltadt*( A*(x[2,:] - x[1,:]) + C*(x[3,:] - x[1,:])) + ...
            #             -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[2])) - 0.5*(this_hsq[1] + this_hsq[2]) / a^4)^2 *(x[2,:] - x[1,:]) + ...
            #             -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[3])) - 0.5*(this_hsq[1] + this_hsq[3]) / c^4)^2 *(x[3,:] - x[1,:])

            t1 = 0.4 / Cdelta^2 * 2 * dCdeltadt * ( A*(x[2,:] - x[1,:]) .+ C*(x[3,:] - x[1,:]))
            t2 = -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (x[2,:] - x[1,:])
            t3 = -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (x[3,:] - x[1,:])

            gradF[:,i] = gradF[:,i] + t1[1] + t2[1] + t3[1]



        end

        tang_proj = I - normals[:,i] * normals[:,i]'
        gradF[:,i] = tang_proj * gradF[:,i]

    end

    return gradF

end

function make_gradF(normals,triangles, vertices, V)

        # normals = normals'
        # triangles = triangles'
        # vertices = vertices'
        # V = V'

        Nvertices = size(vertices, 2)
        gradF = zeros(3,Nvertices)

        for i = 1:Nvertices

            # finds the i-th triangle indices in the triangle 2D array
            mask = findall(x-> x == i, faces)
            num = length(mask)
            (row, col) = zeros(Int64, 1, num), zeros(Int64, 1, num)
            for n in 1:num
                row[n], col[n] = mask[n][1], mask[n][2]
            end


            for this_triang = 1:length(row)
                j1 = (row[this_triang])%(3) + 1 # from 1 to 3
                j2 = (row[this_triang] + 1)%(3) +1 # from 1 to 3


                # radiusvektori uz sho trijsturu virsotneem
                # pirmais atbilst i

                # println()
                # println("triang $i, kaims = $this_triang ----------------------------------------------------------------------------")
                #
                # println()
                # println("x index")
                # println("trig size = ", size(triangles))
                # println("vert size = ", size(vertices))
                #
                # println("$(row[this_triang]), $(col[this_triang]) -> ", triangles[row[this_triang], col[this_triang]])
                # println("$j1, $(col[this_triang]) -> ", triangles[j1, col[this_triang]])
                # println("$j2, $(col[this_triang]) -> ", triangles[j2, col[this_triang]])
                #
                # println(size(triangles))
                # #println(triangles'[1:20, :])
                # println(triangles[1, 1:15])
                # println(triangles[2, 1:15])
                # println(triangles[3, 1:15])
                # println()


                x = [vertices[:,triangles[row[this_triang],col[this_triang]]],
                    vertices[:,triangles[j1, col[this_triang]]],
                    vertices[:,triangles[j2, col[this_triang]]]]

                v = [V[:, triangles[row[this_triang],col[this_triang]]],
                    V[:, triangles[j1, col[this_triang]]],
                    V[:, triangles[j2, col[this_triang]]]]

                # this_hsq = [hsq[triangles[row[this_triang],col[this_triang]]]
                #     hsq[triangles[row[this_triang],j1]]
                #     hsq[triangles[row[this_triang],j2]]]

                a = norm(x[2] - x[1])
                b = norm(x[3] - x[2])
                c = norm(x[1] - x[3])

                Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

                A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)
                B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)
                C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                    (4*Cdelta * (a^2 + b^2 + c^2)^3)

                dCdeltadt = -A * dot(x[2] - x[1], v[2] - v[1]) +
                            -B * dot(x[3] - x[2], v[3] - v[2]) +
                            -C * dot(x[1] - x[3], v[1] - v[3])

                # gradF[i,:] = gradF[i,:] + ...
                #             0.4 / Cdelta^2 * 2 * dCdeltadt*( A*(x[2,:] - x[1,:]) + C*(x[3,:] - x[1,:])) + ...
                #             -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[2])) - 0.5*(this_hsq[1] + this_hsq[2]) / a^4)^2 *(x[2,:] - x[1,:]) + ...
                #             -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[3])) - 0.5*(this_hsq[1] + this_hsq[3]) / c^4)^2 *(x[3,:] - x[1,:])

                t1 = 0.4 / Cdelta^2 * 2 * dCdeltadt * ( A*(x[2,:] - x[1,:]) .+ C*(x[3,:] - x[1,:]))
                t2 = -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (x[2,:] - x[1,:])
                t3 = -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (x[3,:] - x[1,:])

                gradF[:,i] = gradF[:,i] + t1[1] + t2[1] + t3[1]

            end

        end

        return gradF

end

function make_Vvecs_conjgrad(normals,triangles, vertices, vvecs, epsilon, maxIters)

    # [k1, k2] = principal_curvatures[CDE]; # k1 >= k2
    # LAMBDA = k1.^2 + k2.^2 + 0.004
    # K = 4/(sqrt(3) * size(triangles,1)) * sum(LAMBDA.^0.25 .* deltaS)
    # hsq = K * LAMBDA.^(-0.25)
    # triangles = triangles'
    # vertices = vertices'
    # normals = normals'
    # vvecs = vvecs'
    println("passive stabbing")
    # first gradient descent
    f = make_tanggradF(normals,triangles, vertices, vvecs)
    gradFv = make_gradF(normals, triangles, vertices, vvecs)
    gradFf = make_gradF(normals, triangles, vertices, f)

    ksi = - sum(sum(gradFv .* f, dims=2)) / sum(sum(gradFf .* f, dims=2))

    V = vvecs + ksi*f
    Vp = vvecs

    F = make_F(triangles, vertices, V)

    # then conjugated gradient()
    for i = 1:maxIters
        f = make_tanggradF(normals,triangles, vertices, V)
        gradFv = make_gradF(normals, triangles, vertices, V)
        gradFvp = make_gradF(normals, triangles, vertices, Vp)
        gradFf = make_gradF(normals, triangles, vertices, f)
        Ff = make_F(triangles, vertices, f)
        Fdeltav = make_F(triangles, vertices, V-Vp)

        a1 = sum(sum(gradFv .* f, dims=2))
        b1 = Ff
        c1 = sum(sum((gradFv .- gradFvp) .* f, dims=2))
        a2 = sum(sum(gradFv .* (V.-Vp), dims=2))
        b2 = sum(sum(gradFf .* (V.-Vp), dims=2))
        c2 = Fdeltav

        ksi = (a2*c1 - 2*a1*c2) / (4*b1*c2 - b2*c1)
        eta = (2*a2*b1 - a1*b2) / (b2*c1 - 4*b1*c2)

        Vtemp = V
        V = V .+ ksi*f .+ eta*(V.-Vp)
        Vp = Vtemp

        Fp = F
        F = make_F(triangles, vertices, V)
        #println(F)
        #println((Fp-F)/F)

        if (Fp - F)/F < epsilon
            println("improved tangential velocities")
            break
        end
        if i == maxIters
            println("tangential velocities not fully converged")
        end
    end

    return V
end

function passive_stab(normals,triangles, vertices, vvecs, epsilon, maxIters)
    # [k1, k2] = principal_curvatures[CDE]; # k1 >= k2
    # LAMBDA = k1.^2 + k2.^2 + 0.004
    # K = 4/(sqrt(3) * size(triangles,1)) * sum(LAMBDA.^0.25 .* deltaS)
    # hsq = K * LAMBDA.^(-0.25)
    # triangles = triangles'
    # vertices = vertices'
    # normals = normals'
    # vvecs = vvecs'

    println("passive stabbing")
    # first gradient descent
    f = make_tanggradF(normals,triangles, vertices, vvecs)
    gradFv = make_gradF(normals, triangles, vertices, vvecs)
    gradFf = make_gradF(normals, triangles, vertices, f)

    ksi = - sum(sum(gradFv .* f, dims=2)) / sum(sum(gradFf .* f, dims=2))

    V = vvecs + ksi*f
    Vp = vvecs

    F = make_F(triangles, vertices, V)

    # then conjugated gradient()
    for i = 1:maxIters
        f = make_tanggradF(normals,triangles, vertices, V)
        gradFv = make_gradF(normals, triangles, vertices, V)
        gradFvp = make_gradF(normals, triangles, vertices, Vp)
        gradFf = make_gradF(normals, triangles, vertices, f)
        Ff = make_F(triangles, vertices, f)
        Fdeltav = make_F(triangles, vertices, V-Vp)

        a1 = sum(sum(gradFv .* f, dims=2))
        b1 = Ff
        c1 = sum(sum((gradFv .- gradFvp) .* f, dims=2))
        a2 = sum(sum(gradFv .* (V.-Vp), dims=2))
        b2 = sum(sum(gradFf .* (V.-Vp), dims=2))
        c2 = Fdeltav

        ksi = (a2*c1 - 2*a1*c2) / (4*b1*c2 - b2*c1)
        eta = (2*a2*b1 - a1*b2) / (b2*c1 - 4*b1*c2)

        Vtemp = V
        V = V .+ ksi*f .+ eta*(V.-Vp)
        Vp = Vtemp

        Fp = F
        F = make_F(triangles, vertices, V)

        if (Fp - F)/F < epsilon
            println("improved tangential velocities")
            break
        end
        if i == maxIters
            println("tangential velocities not fully converged")
        end
    end

    return V
end

function make_enright_velocities(points, t)
    x, y, z = points[1, :], points[2, :], points[3, :]

    vx = 2*sin.(pi*x).^2 .* sin.(2pi*y) .* sin.(2pi*z) .* sin(2/3*pi*t)
    vy = -sin.(2pi*x) .* sin.(pi*y).^2 .* sin.(2pi*z) .* sin(2/3*pi*t)
    vz = -sin.(2pi*x) .* sin.(2pi*y) .* sin.(pi*z).^2 .* sin(2/3*pi*t)

    v = [vx vy vz]

    return Array{Float64}(transpose(v))
end

function NeighborVertices(vertex, faces)
    mask = findall(x-> x == vertex, faces)
    neighbors = zeros(Int64, 3, length(mask))
    num_neighbors = length(mask)
    for n in 1:num_neighbors
        neighbors[:, n] = faces[:,mask[n][2]]
    end
    return setdiff(neighbors, vertex)
end

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
        # println("checking if many neighbors are marked")
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
                    marked_faces[i] = true
                    # println(marked_faces[i])
                    # println(nearby_marked_faces," marked kaimiņi")
                end

            end
        end

        if marked_faces_old == marked_faces
            neighbors_marked = true
        else
            println("marked new faces")
        end
    end # end while


    return marked_faces
end

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

function active_stabilize_old_surface(points_old,CDE_old,normals_old,points0::Array{Float64,2},faces::Array{Int64,2},connectivity::Array{Int64,2},edges::Array{Int64,2};
    deltakoef=0.01, R0=1.0, gamma=0.25, p=50, r=100, checkiters=100, maxiters=1000,critSc = 0.75,critCdelta = 1.15)
    # actively rearange vertices on a surfaces given by fitted paraboloids
    # this has been modified to use points before splitting for surface determination
    # as per Zinchenko(2013)
    println("active stabilization")
    points = copy(points0)

    closefaces = make_closefaces(faces)
    dS = make_dS(points,faces)

    k1 = Array{Float64}(undef, size(points,2))
    k2 = Array{Float64}(undef, size(points,2))
    for i = 1:size(points,2)
        r0 = points[:,i]
        minind = argmin(sum((points_old .- r0).^2,dims=1))[2]
        r0 = to_local(r0-points_old[:,minind],normals_old[:,minind])
        k1[i],k2[i] =  make_pc_local(CDE_old[:,minind],r0[1],r0[2])
    end

    LAMBDA = k1.^2 + k2.^2 .+ 0.004/R0^2
    K = 4/(sqrt(3) * size(faces,2)) * sum(LAMBDA.^gamma .* dS)
    hsq = K * LAMBDA.^(-gamma)

    no_improvement = true
    # initiallize improvement criteria
    xij = make_edge_lens(points,edges)
    hij = sqrt.(0.5*(hsq[edges[1,:]].^2 + hsq[edges[2,:]].^2))

    Sc0 = maximum(xij./hij) / minimum(xij./hij)
    Cdelta_min0 = minimum(make_Cdeltas(points, faces))
    for iter = 1:maxiters
        #println(iter)
        gradE = make_gradE(points,faces,closefaces,hsq; p=p,r=r)
        delta = deltakoef * minimum(make_min_edges(points,connectivity) ./ sum(sqrt.(gradE.^2),dims=1))

        #println("dPoints= ",delta*maximum(sum(sqrt.(gradE.^2),dims=1)))
        #println("E = ", make_E(points,faces,hsq; p=p,r=r))


        points = points - delta * gradE

        #project points on the drop
        for i = 1:size(points,2)
            points[:,i] = project_on_drop(points_old,CDE_old,normals_old,points[:,i])
        end

        # recalculate the mesh parameters
        dS = make_dS(points,faces)
        #k1,k2 = make_pc(CDE)
        for i = 1:size(points,2)
            r0 = points[:,i]
            minind = argmin(sum((points_old .- r0).^2,dims=1))[2]
            r0 = to_local(r0-points_old[:,minind],normals_old[:,minind])
            k1[i],k2[i] =  make_pc_local(CDE_old[:,minind],r0[1],r0[2])
        end
        LAMBDA = k1.^2 + k2.^2 .+ 0.004/R0^2
        K = 4/(sqrt(3) * size(faces,2)) * sum(LAMBDA.^gamma .* dS)
        hsq = K * LAMBDA.^(-gamma)
        if iter < checkiters
            xij = make_edge_lens(points,edges)
            hij = sqrt.(0.5*(hsq[edges[1,:]].^2 + hsq[edges[2,:]].^2))

            Sc = maximum(xij./hij) / minimum(xij./hij)
            Cdelta_min = minimum(make_Cdeltas(points, faces))

            #println("Sc/Sc0 = ",Sc/Sc0)
            #println("Cdelta/Cdelta0 = ",Cdelta_min/Cdelta_min0)
            if Sc > critSc*Sc0 && Cdelta_min < critCdelta*Cdelta_min0
                no_improvement = false
            end
        end

        if iter == checkiters
            if no_improvement == true
                println("no significant improvement achieved")
                println("reversing changes")
                points = points0
                break
            else
                println("improvement detected in the first ", checkiters, " iterations")
                println("iterating for ", maxiters - checkiters, " more iterations")
            end
        end

        # e = make_E(points, faces, hsq,p=p, r=r)
        # println("E = $e")

        if iter%500 == 0
            println("iteration ",iter)

        end
    end
    return points
end
