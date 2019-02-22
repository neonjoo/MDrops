function make_magvelocities(vertices, normals, lambda, Bm, mu, Hn_2, Ht_2)
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
