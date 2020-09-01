function make_magvelocities_old(vertices, normals, lambda, Bm, mu, Hn_2, Ht_2)
    # lambda = int viscosity / ext viscosity
    # Returns vertex velocities from the Stokes flow integral equations
    # using Wielandt's deflation as describeed in Zinchenko 1997, Das 2017
    # w + T*w + 4*pi*C*w + B*w = F
    println("calculating velocities")
    deltaS = make_dS(points,faces)

    #fbar = mu*(mu-1)/8pi * Hfieldn.^2 + (mu-1)/8pi * Hfieldt.^2
    fbar = mu*(mu-1)/8pi * Hn_2 + (mu-1)/8pi * Ht_2

    S = sum(deltaS)
    xc = [0.,0.,0.]
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        xc = xc + rx * deltaS[x] / S
    end

    F = zeros(3*size(vertices,2),1)
    T = zeros(3*size(vertices,2),3*size(vertices,2))
    #C = C1+C2
    C1 = zeros(3*size(vertices,2),3*size(vertices,2))
    C2 = zeros(3*size(vertices,2),3*size(vertices,2))
    B = zeros(3*size(vertices,2),3*size(vertices,2))

    M = zeros(3,3)
    # fill M
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        # Matrix(1.0I, 3, 3) <- 3x3 Identity matrix
        M += (Matrix(1.0I, 3, 3) * dot(rx - xc, rx - xc) -
            (rx - xc) * (rx - xc)') * deltaS[x]
            #   ^ changed the transpose
    end
    Minv = inv(M) # calculate inverse of M

    #create Levi-Civita symbol
    epsilon = cat(dims=3, [0 0  0;  0 0 1; 0 -1 0],
                              [0 0 -1;  0 0 0; 1  0 0],
                              [0 1  0; -1 0 0; 0  0 0])


    for y = 1:size(vertices,2)
        ry = vertices[:,y]

        for j = 1:3

            for x = 1:size(vertices,2)
                rx = vertices[:,x]

                Gji = (Matrix(1.0I, (3,3)) / norm(rx-ry) +
                (rx-ry) * (rx-ry)' / norm(rx-ry)^3)

                if x != y
                   F[3*y-3+j] +=
                       1/(8*pi) *(
                       dot(rx-ry,normals[:,x]) * normals[j,y]+
                       dot(rx-ry,normals[:,y]) * normals[j,x]+
                       (1-dot(normals[:,x],normals[:,y])) * (rx[j]-ry[j])-
                       3*(rx[j]-ry[j])/norm(rx-ry)^2*
                       dot(normals[:,x]+normals[:,y],rx-ry)*
                       dot(rx-ry,normals[:,x])) * deltaS[x]/norm(rx-ry)^3

                # magnetic part
                   F[3*y-3+j] += Bm/(8*pi) *
                   (fbar[x] - fbar[y]) *
                   dot(Gji[j,:], normals[:,x]) * deltaS[x]
                #else
                #polar integration scheme
                end

                for i = 1:3
                    if x != y
                       T[3*y-3+j,3*x-3+i] = -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]

                       T[3*y-3+j,3*y-3+i] -= -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]
                    end

                    B[3*y-3+j,3*x-3+i] = -4*pi/S * normals[j,y] * normals[i,x] * deltaS[x]

                    C1[3*y-3+j,3*x-3+i] =  deltaS[x]/S

                    for k = 1:3
                    for m = 1:3
                    for n = 1:3
                    for p = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            epsilon[j,k,m]*Minv[k,n]*epsilon[n,p,i]*
                            (rx[p] - xc[p])*deltaS[x]*(ry[m] - xc[m])
                    end
                    end
                    end
                    end

                end
            end

        end
    end
    C = C1 + C2

    # w + T*w + 4*pi*C*w + B*w = F
    #w = (eye(size(T)) + T + 4*pi*C + B) \ F; old and bad
    w = (Matrix(1.0I, size(T)) + T + (lambda-1)/2*C + (lambda-1)/(8*pi)*B) \ F

    # w' = C*w
    # v = w + (lambda-1)/2 * w'
    varr = w + (lambda-1)/2 * C*w


    velocities = reshape(varr ,3,size(vertices,2))
    return velocities
end

function make_magvelocities(vertices, normals, lambda, Bm, mu, Hn_2, Ht_2)
    #without the use of levi civita tensor
    # lambda = int viscosity / ext viscosity
    # Returns vertex velocities from the Stokes flow integral equations
    # using Wielandt's deflation as describeed in Zinchenko 1997, Das 2017
    # w + T*w + 4*pi*C*w + B*w = F
    println("calculating velocities")
    deltaS = make_dS(points,faces)

    #fbar = mu*(mu-1)/8pi * Hfieldn.^2 + (mu-1)/8pi * Hfieldt.^2
    fbar = mu*(mu-1)/8pi * Hn_2 + (mu-1)/8pi * Ht_2


    S = sum(deltaS)
    xc = [0.,0.,0.]
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        xc = xc + rx * deltaS[x] / S
    end

    F = zeros(3*size(vertices,2),1)
    T = zeros(3*size(vertices,2),3*size(vertices,2))
    #C = C1+C2
    C1 = zeros(3*size(vertices,2),3*size(vertices,2))
    C2 = zeros(3*size(vertices,2),3*size(vertices,2))
    B = zeros(3*size(vertices,2),3*size(vertices,2))

    M = zeros(3,3)
    # fill M
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        # Matrix(1.0I, 3, 3) <- 3x3 Identity matrix
        M += (Matrix(1.0I, 3, 3) * dot(rx - xc, rx - xc) -
            (rx - xc) * (rx - xc)') * deltaS[x]
            #   ^ changed the transpose
    end
    Minv = inv(M) # calculate inverse of M

    for y = 1:size(vertices,2)
        ry = vertices[:,y]

        for j = 1:3

            for x = 1:size(vertices,2)
                rx = vertices[:,x]

                Gji = (Matrix(1.0I, (3,3)) / norm(rx-ry) +
                (rx-ry) * (rx-ry)' / norm(rx-ry)^3)

                if x != y
                   F[3*y-3+j] +=
                       1/(8*pi) *(
                       dot(rx-ry,normals[:,x]) * normals[j,y]+
                       dot(rx-ry,normals[:,y]) * normals[j,x]+
                       (1-dot(normals[:,x],normals[:,y])) * (rx[j]-ry[j])-
                       3*(rx[j]-ry[j])/norm(rx-ry)^2*
                       dot(normals[:,x]+normals[:,y],rx-ry)*
                       dot(rx-ry,normals[:,x])) * deltaS[x]/norm(rx-ry)^3

                # magnetic part
                   F[3*y-3+j] += Bm/(8*pi) *
                   (fbar[x] - fbar[y]) *
                   dot(Gji[j,:], normals[:,x]) * deltaS[x]
                #else
                #polar integration scheme
                end

                for i = 1:3
                    if x != y
                       T[3*y-3+j,3*x-3+i] = -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]

                       T[3*y-3+j,3*y-3+i] -= -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]
                    end

                    B[3*y-3+j,3*x-3+i] = -4*pi/S * normals[j,y] * normals[i,x] * deltaS[x]

                    C1[3*y-3+j,3*x-3+i] =  deltaS[x]/S

                    for k = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            Minv[k,j]*(rx[k] - xc[k])*deltaS[x]*(ry[i] - xc[i]) -
                            Minv[i,j]*(rx[k] - xc[k])*deltaS[x]*(ry[k] - xc[k]) -
                            Minv[k,k]*(rx[j] - xc[j])*deltaS[x]*(ry[i] - xc[i]) +
                            Minv[i,k]*(rx[j] - xc[j])*deltaS[x]*(ry[k] - xc[k])
                    end
                    if i==j
                        for k = 1:3
                        for m = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            Minv[k,k]*(rx[m] - xc[m])*deltaS[x]*(ry[m] - xc[m]) -
                            Minv[k,m]*(rx[k] - xc[k])*deltaS[x]*(ry[m] - xc[m])
                        end
                        end
                    end

                end
            end

        end
    end
    C = C1 + C2

    # w + T*w + 4*pi*C*w + B*w = F
    #w = (eye(size(T)) + T + 4*pi*C + B) \ F; old and bad
    w = (Matrix(1.0I, size(T)) + T + (lambda-1)/2*C + (lambda-1)/(8*pi)*B) \ F

    # w' = C*w
    # v = w + (lambda-1)/2 * w'
    varr = w + (lambda-1)/2 * C*w


    velocities = reshape(varr ,3,size(vertices,2))
    return velocities
end


function make_magvelocities_no_Wie(vertices, normals, lambda, Bm, mu, Hn_2, Ht_2)
    # Returns vertex velocities from the Stokes flow integral equations without Wielandt's deflation
    # v = Av + F
    println("calculation velocities, no Wieldandt")
    deltaS = make_dS(points,faces)

    #fbar = mu*(mu-1)/8pi * Hfieldn.^2 + (mu-1)/8pi * Hfieldt.^2
    fbar = mu*(mu-1)/8pi * Hn_2 + (mu-1)/8pi * Ht_2

    F = zeros(3*size(vertices,2),1)
    A = zeros(3*size(vertices,2), 3*size(vertices,2))

    for y = 1:size(vertices,2)
        ry = vertices[:,y]

        for j = 1:3

            for x = 1:size(vertices,2)
                rx = vertices[:,x]

                Gji = (Matrix(1.0I, (3,3)) / norm(rx-ry) +
                (rx-ry) * (rx-ry)' / norm(rx-ry)^3)

                if x != y
                   F[3*y-3+j] +=
                       1/(8*pi) *(
                       dot(rx-ry,normals[:,x]) * normals[j,y]+
                       dot(rx-ry,normals[:,y]) * normals[j,x]+
                       (1-dot(normals[:,x],normals[:,y])) * (rx[j]-ry[j])-
                       3*(rx[j]-ry[j])/norm(rx-ry)^2*
                       dot(normals[:,x]+normals[:,y],rx-ry)*
                       dot(rx-ry,normals[:,x])) * deltaS[x]/norm(rx-ry)^3

                # magnetic part
                   F[3*y-3+j] += Bm/(8*pi) *
                   (fbar[x] - fbar[y]) *
                   dot(Gji[j,:], normals[:,x]) * deltaS[x]
                end

                # build the A matrix
                for i = 1:3
                    if x != y
                        A[3y-3+j, 3x-3+i] = (1-lambda)/8/pi * (-6 * (rx[i]-ry[i])*(rx[j]-ry[j])) / norm(rx-ry)^5 *
                            dot(rx-ry, normals[:,x]) * deltaS[x]
                        A[3y-3+j, 3y-3+i] += (1-lambda)/8/pi * (-6 * (rx[i]-ry[i])*(rx[j]-ry[j])) / norm(rx-ry)^5 *
                            dot(rx-ry, normals[:,x]) * deltaS[x]


                    end # if x!=y
                end # i for
            end # x for
        end # j for
    end # y for

    varr = (A-I) \ F;  # solve integral eq. for velocity

    velocities = reshape(varr ,3,size(vertices,2))
    return velocities

end # function

function make_magvelocities_2(vertices, normals, lambda, Bm, mu, Hn_2, Ht_2)
    # early exit for lambda=1
    #without the use of levi civita tensor
    # lambda = int viscosity / ext viscosity
    # Returns vertex velocities from the Stokes flow integral equations
    # using Wielandt's deflation as describeed in Zinchenko 1997, Das 2017
    # w + T*w + 4*pi*C*w + B*w = F
    println("calculating velocities")
    deltaS = make_dS(points,faces)

    #fbar = mu*(mu-1)/8pi * Hfieldn.^2 + (mu-1)/8pi * Hfieldt.^2
    fbar = mu*(mu-1)/8pi * Hn_2 + (mu-1)/8pi * Ht_2


    S = sum(deltaS)
    xc = [0.,0.,0.]
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        xc = xc + rx * deltaS[x] / S
    end

    F = zeros(3*size(vertices,2),1)
    T = zeros(3*size(vertices,2),3*size(vertices,2))
    #C = C1+C2
    C1 = zeros(3*size(vertices,2),3*size(vertices,2))
    C2 = zeros(3*size(vertices,2),3*size(vertices,2))
    B = zeros(3*size(vertices,2),3*size(vertices,2))



    if lambda == 1
        for y = 1:size(vertices,2)
            ry = vertices[:,y]

            for j = 1:3

                for x = 1:size(vertices,2)
                    rx = vertices[:,x]

                    Gji = (Matrix(1.0I, (3,3)) / norm(rx-ry) +
                    (rx-ry) * (rx-ry)' / norm(rx-ry)^3)

                    if x != y
                       F[3*y-3+j] +=
                           1/(8*pi) *(
                           dot(rx-ry,normals[:,x]) * normals[j,y]+
                           dot(rx-ry,normals[:,y]) * normals[j,x]+
                           (1-dot(normals[:,x],normals[:,y])) * (rx[j]-ry[j])-
                           3*(rx[j]-ry[j])/norm(rx-ry)^2*
                           dot(normals[:,x]+normals[:,y],rx-ry)*
                           dot(rx-ry,normals[:,x])) * deltaS[x]/norm(rx-ry)^3

                    # magnetic part
                       F[3*y-3+j] += Bm/(8*pi) *
                       (fbar[x] - fbar[y]) *
                       dot(Gji[j,:], normals[:,x]) * deltaS[x]
                    #else
                    #polar integration scheme
                    end
                end

            end
        end
        # w + T*w + 4*pi*C*w + B*w = F
        #w = (eye(size(T)) + T + 4*pi*C + B) \ F; old and bad
        w = (Matrix(1.0I, size(T)) + T) \ F

        # w' = C*w
        # v = w + (lambda-1)/2 * w'
        varr = w


        velocities = reshape(varr ,3,size(vertices,2))
        return velocities
    end



    M = zeros(3,3)
    # fill M
    for x = 1:size(vertices,2)
        rx = vertices[:,x]
        # Matrix(1.0I, 3, 3) <- 3x3 Identity matrix
        M += (Matrix(1.0I, 3, 3) * dot(rx - xc, rx - xc) -
            (rx - xc) * (rx - xc)') * deltaS[x]
            #   ^ changed the transpose
    end
    Minv = inv(M) # calculate inverse of M

    for y = 1:size(vertices,2)
        ry = vertices[:,y]

        for j = 1:3

            for x = 1:size(vertices,2)
                rx = vertices[:,x]

                Gji = (Matrix(1.0I, (3,3)) / norm(rx-ry) +
                (rx-ry) * (rx-ry)' / norm(rx-ry)^3)

                if x != y
                   F[3*y-3+j] +=
                       1/(8*pi) *(
                       dot(rx-ry,normals[:,x]) * normals[j,y]+
                       dot(rx-ry,normals[:,y]) * normals[j,x]+
                       (1-dot(normals[:,x],normals[:,y])) * (rx[j]-ry[j])-
                       3*(rx[j]-ry[j])/norm(rx-ry)^2*
                       dot(normals[:,x]+normals[:,y],rx-ry)*
                       dot(rx-ry,normals[:,x])) * deltaS[x]/norm(rx-ry)^3

                # magnetic part
                   F[3*y-3+j] += Bm/(8*pi) *
                   (fbar[x] - fbar[y]) *
                   dot(Gji[j,:], normals[:,x]) * deltaS[x]
                #else
                #polar integration scheme
                end

                for i = 1:3
                    if x != y
                       T[3*y-3+j,3*x-3+i] = -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]

                       T[3*y-3+j,3*y-3+i] -= -(1-lambda)/(8*pi)*
                           (-6*(rx[i]-ry[i])*(rx[j]-ry[j]) / norm(rx-ry)^5 )*
                           dot(rx-ry,normals[:,x])*
                           deltaS[x]
                    end

                    B[3*y-3+j,3*x-3+i] = -4*pi/S * normals[j,y] * normals[i,x] * deltaS[x]

                    C1[3*y-3+j,3*x-3+i] =  deltaS[x]/S

                    for k = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            Minv[k,j]*(rx[k] - xc[k])*deltaS[x]*(ry[i] - xc[i]) -
                            Minv[i,j]*(rx[k] - xc[k])*deltaS[x]*(ry[k] - xc[k]) -
                            Minv[k,k]*(rx[j] - xc[j])*deltaS[x]*(ry[i] - xc[i]) +
                            Minv[i,k]*(rx[j] - xc[j])*deltaS[x]*(ry[k] - xc[k])
                    end
                    if i==j
                        for k = 1:3
                        for m = 1:3
                        C2[3*y-3+j,3*x-3+i] +=
                            Minv[k,k]*(rx[m] - xc[m])*deltaS[x]*(ry[m] - xc[m]) -
                            Minv[k,m]*(rx[k] - xc[k])*deltaS[x]*(ry[m] - xc[m])
                        end
                        end
                    end

                end
            end

        end
    end
    C = C1 + C2

    # w + T*w + 4*pi*C*w + B*w = F
    #w = (eye(size(T)) + T + 4*pi*C + B) \ F; old and bad
    w = (Matrix(1.0I, size(T)) + T + (lambda-1)/2*C + (lambda-1)/(8*pi)*B) \ F

    # w' = C*w
    # v = w + (lambda-1)/2 * w'
    varr = w + (lambda-1)/2 * C*w


    velocities = reshape(varr ,3,size(vertices,2))
    return velocities
end

function scalar_magnetic_potential(r::Array{Float64}, M::Float64, angle::Float64,  xlims::Array{Float64}, ylims::Array{Float64}, zlims::Array{Float64})
    # Scalar magnetic potential of a block magnet, bounded by lims (according to https://contattafiles.s3.us-west-1.amazonaws.com/tnt41611/uDDJpl0uRpMsYLf/Ravaud_Lemarquand_2009_Magnetic%20field%20produced%20by%20a%20parallelepipedic%20magnet%20of%20various%20and%20uniform.pdf)
    # R. Ravaud and G. Lemarquand "MAGNETIC FIELD PRODUCED BY A PARALLELEPIPEDIC MAGNET OF VARIOUS AND UNIFORM POLARIZATION"
    # M - magnetization, angle wrt to X axis in XY plane

    function dist(i::Int64, j::Int64, k::Int64, r::Array{Float64}, xlims::Array{Float64}, ylims::Array{Float64}, zlims::Array{Float64})
        # distance function wrt to given faces of the block magnet
        x, y, z = r
        xi, yj, zk = xlims[i], ylims[j], zlims[k]

        return sqrt((x-xi)^2 + (y-yj)^2 + (z-zk)^2)
    end

    function phi_1(i::Int64, j::Int64, k::Int64, r::Array{Float64}, xlims::Array{Float64}, ylims::Array{Float64}, zlims::Array{Float64})
        # helper function, see ref above
        x, y, z = r
        xi, yj, zk = xlims[i], ylims[j], zlims[k]
        d = dist(i, j, k, r, xlims, ylims, zlims)

        return zk + (x-xi)*atan((z-zk)/(x-xi)) + (y-yj)*log(z-zk+d) -
                    (x-xi)*atan( (y-yj)*(z-zk) / ((x-xi)*d) ) + (z-zk)*log(y-yj+d)
    end

    function phi_2(i::Int64, j::Int64, k::Int64, r::Array{Float64}, xlims::Array{Float64}, ylims::Array{Float64}, zlims::Array{Float64})
        # helper function, see ref above
        x, y, z = r
        xi, yj, zk = xlims[i], ylims[j], zlims[k]
        d = dist(i, j, k, r, xlims, ylims, zlims)

        return zk + (y-yj)*atan((z-zk)/(y-yj)) + (x-xi)*log(z-zk+d) -
                    (y-yj)*atan( (x-xi)*(z-zk) / ((y-yj)*d) ) + (z-zk)*log(x-xi+d)
    end

    pot = 0
    for i in 1:2
        for j in 1:2
            for k in 1:2
                pot = pot + (-1)^(i+j+k) * ( cos(angle)*phi_1(i, j, k, r, xlims, ylims, zlims) + sin(angle)*phi_2(i, j, k, r, xlims, ylims, zlims) )

            end
        end
    end

    return pot*M/4pi/(4pi*10^-7)

end

function quadropole_potential(point::Array{Float64})
    x, y, z = point
    M = 1.0
    z_lim = 1.75
    return scalar_magnetic_potential([x, y, z], M, 1pi, [-1.0, -0.5], [0., 0.5], [-z_lim/2, z_lim/2]) +
                 scalar_magnetic_potential([x, y, z], M, 0., [-1.0, -0.5], [-0.5, 0.], [-z_lim/2, z_lim/2]) +
                 scalar_magnetic_potential([x, y, z], M, 1pi, [0.5, 1.], [0., 0.5], [-z_lim/2, z_lim/2]) +
                 scalar_magnetic_potential([x, y, z], M, 0., [0.5, 1.], [-0.5, 0.], [-z_lim/2, z_lim/2])
end

function make_L_sing(points, faces, normals; gaussorder=5)
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = zeros(Float64, N)

    function make_n(x,x1,x2,x3,n1,n2,n3) # normal linear interpolation
        A = [x1 x2 x3] # matrix of vertex radiusvecotrs
        B = [n1 n2 n3] # matrix of vertex normals

        zeta_xi_eta = A \ x # find local triangle parameters

        n = B * zeta_xi_eta
        return n/norm(n)
    end


    for ykey in 1:N
        function make_L(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey],ny=normals[:,ykey])
            nx = make_n(x,x1,x2,x3,n1,n2,n3)
            r = y - x

            return dot(r, nx-ny) / norm(r)^3 / 4pi
        end

        function make_L_times_norm_r(x,x1,x2,x3,n1,n2,n3; y=points[:,ykey], ny=normals[:,ykey])
            nx = make_n(x,x1,x2,x3,n1,n2,n3)
            r = y - x

            return dot(r, nx-ny) / norm(r)^2 / 4pi
        end

        ny = normals[:, ykey]
        ry = points[:, ykey]

        for i in 1:size(faces, 2)
            #println("triangle: $i")
            if !(ykey in faces[:,i]) # if not singular triangle
                x1, x2, x3 = [points[:, faces[j,i]] for j in 1:3]
                n1, n2, n3 = [normals[:, faces[j,i]] for j in 1:3]

                L[ykey] += gauss_nonsingular(x->make_L(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #println(x1, x2, x3)
            else # if is singular triangle
                singul_ind = findfirst(ind->ind==ykey,faces[:,i])

                x1 = points[:,faces[singul_ind,i]]
                x2 = points[:,faces[(singul_ind) % 3 + 1,i]]
                x3 = points[:,faces[(singul_ind + 1) % 3 + 1,i]]

                n1 = normals[:,faces[singul_ind,i]]
                n2 = normals[:,faces[(singul_ind) % 3 + 1,i]]
                n3 = normals[:,faces[(singul_ind + 1) % 3 + 1,i]]
                #println("L = $L,eeee $(L[ykey])")
                #L[ykey] +=
                L[ykey] += gauss_weaksingular(x->make_L_times_norm_r(x,x1,x2,x3,n1,n2,n3), x1,x2,x3,gaussorder)
                #println("$Ltemp, $(L[ykey]), $ykey")
            end

        end
    end

    return -L # implement the minus in the end. For a sphere L=-1 for all points
end

function make_deltaH_normal(points, faces, normals, mu, H0; gaussorder=5)
    # calculate the jump in the normal H field on drop surface
    # from D. Das and D. Saintillan "Electrohydrodynamics of viscous drops in strong electric fields: numerical simulations"
    alpha = mu
    deltaS = make_dS(points,faces)
    N = size(points, 2)
    L = make_L_sing(points, faces, normals; gaussorder=gaussorder)
    L = Diagonal(alpha/(alpha-1) .- L) # alpha term should be positive
    F = zeros(Float64, N, N) # integral equation matrix
    H = zeros(Float64, N)

    for ykey in 1:N
        ny = normals[:, ykey]
        ry = points[:, ykey]

        for xkey in 1:N
            if xkey == ykey
                continue
            end
            nx = normals[:, xkey]
            rx = points[:, xkey]
            F[ykey, xkey] +=  dot(ny, ry-rx) / norm(ry-rx)^3 * deltaS[xkey] / (4pi)

        end
    end

    for ykey in 1:N
        H[ykey] = dot(H0, normals[:, ykey]) # no minus
    end

    F_reg = F - Diagonal([sum(F[i, :]) for i in 1:N])
    deltaH_normal = (L - F_reg) \ H

    return deltaH_normal'
end
