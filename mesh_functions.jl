function make_pc(CDE::Array{Float64,2})
# returns principal curvatures k1, k2 at the vertices
# with locally fitted paraboloids z = Cx^2 + Dxy + Ey^2
# k1 > k2

    C = CDE[1,:]
    D = CDE[2,:]
    E = CDE[3,:]

    k1 = -C - E + sqrt.(C.^2 + D.^2 - 2*C.*E + E.^2)
    k2 = -C - E - sqrt.(C.^2 + D.^2 - 2*C.*E + E.^2)

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
        0.5*sqrt(
            4*(D^2 - 4*C*E)*(1 + 4*C^2*x^2 + 4*C*D*x*y + 4*D*E*x*y + 4*E^2*y^2 + D^2*(x^2 + y^2))^2 +
            4*(1 + (2*C*x + D*y)^2 + (D*x + 2*E*y)^2) * (E + 4*C^2*E*x^2 - D^3*x*y - D^2*E*y^2 + C*(1 - D^2*x^2 + 4*D*E*x*y + 4*E^2*y^2))^2
            )
        )

    k1 = -1/magN^2 *
        (
        C*sqrt(magN) + E*sqrt(magN) - C*D^2*x^2*sqrt(magN) + 4*C^2*E*x^2*sqrt(magN) -
        D^3*x*y*sqrt(magN) + 4*C*D*E*x*y*sqrt(magN) - D^2*E*y^2*sqrt(magN) + 4*C*E^2*y^2*sqrt(magN) -
        0.5*sqrt(
            4*(D^2 - 4*C*E)*(1 + 4*C^2*x^2 + 4*C*D*x*y + 4*D*E*x*y + 4*E^2*y^2 + D^2*(x^2 + y^2))^2 +
            4*(1 + (2*C*x + D*y)^2 + (D*x + 2*E*y)^2) * (E + 4*C^2*E*x^2 - D^3*x*y - D^2*E*y^2 + C*(1 - D^2*x^2 + 4*D*E*x*y + 4*E^2*y^2))^2
            )
        )

    return k1,k2
end

function to_local(r::Array{Float64,1},normal::Array{Float64,1})
    # rotate a vector to local coordinate system
    # with z axis along a normal
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
        if maximum(sqrt.(sum(x -> x^2, normalsp - normals, dims=1))) < eps_outer
            # biggest absolute change in normal vector
            println("paraboloid fit converged")
            break
        end
    end
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
    # ,Optim.Options(f_tol=1.e-10)
    res = Optim.optimize(f,g!,h!,x0,NewtonTrustRegion())
    #println(res)
    x1 = Optim.minimizer(res)[1]
    x2 = Optim.minimizer(res)[2]
    r0 = [x1,x2,CDE[1,i]*x1^2 + CDE[2,i]*x1*x2 + CDE[3,i]*x2^2]
    #println(r0)
    #println(x0)

    r0 = to_global(r0,normals[:,i]) + r
    return r0
end

function active_stabilize(points0::Array{Float64,2},faces::Array{Int64,2},CDE::Array{Float64,2},connectivity::Array{Int64,2},normals::Array{Float64,2};
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
            k1[i],k2[i] =  make_pc_local(CDE[:,i],r0[1],r0[2])
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
            if Sc > critSc*Sc0 || Cdelta_min < critCdelta*Cdelta_min0
                no_improvement = false
            end
        end

        if iter == checkiters
            if no_improvement == true
                println("no significant improvement achieved")
                println("reversing changes")
                #points = points0
                #break
            else
                println("improvement detected in the first ", checkiters, " iterations")
                println("iterating for ", maxiters - checkiters, " more iterations")
            end
        end
        if iter%100 == 0
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
