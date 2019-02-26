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
            break
        end
    end

    return V
end
