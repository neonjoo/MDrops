function make_F(triangles, vertices, V)

    # triangles = triangles'
    # vertices = vertices'
    # V = V'

    Ntriangles = size(triangles,1)

    F = 0

    for i = 1:Ntriangles
        x = [vertices[triangles[i,1],:]
             vertices[triangles[i,2],:]
             vertices[triangles[i,3],:]]

        v = [V[triangles[i,1],:]
             V[triangles[i,2],:]
             V[triangles[i,3],:]]

        # this_hsq = [hsq[triangles[i,1]]
        #             hsq[triangles[i,2]]
        #             hsq[triangles[i,3]]]

        a = norm(x[2,:] - x[1,:])
        b = norm(x[3,:] - x[2,:])
        c = norm(x[1,:] - x[3,:])
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

        dCdeltadt = -A * dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) +
                    -B * dot(x[3,:] - x[2,:], v[3,:] - v[2,:]) +
                    -C * dot(x[1,:] - x[3,:], v[1,:] - v[3,:])

        println(dCdeltadt)

        # F = F + ...
        #     0.4 / Cdelta^2 * dCdeltadt^2 + ...
        #     2*( dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[2])) - 0.5*(this_hsq[1] + this_hsq[2]) / a^4) )^2 + ...
        #     2*( dot(x[3,:] - x[2,:], v[3,:] - v[2,:]) * (1/(0.5*(this_hsq[3] + this_hsq[2])) - 0.5*(this_hsq[3] + this_hsq[2]) / b^4) )^2 + ...
        #     2*( dot(x[1,:] - x[3,:], v[1,:] - v[3,:]) * (1/(0.5*(this_hsq[1] + this_hsq[3])) - 0.5*(this_hsq[1] + this_hsq[3]) / c^4) )^2

        F = F +
            0.4 / Cdelta^2 * dCdeltadt^2 +
            2*( dot(x[2,:] - x[1,:], v[2,:] - v[1,:]))^2 +
            2*( dot(x[3,:] - x[2,:], v[3,:] - v[2,:]))^2 +
            2*( dot(x[1,:] - x[3,:], v[1,:] - v[3,:]))^2

    end

    return F
end



function make_tanggradF(normals,triangles, vertices, V)

    #
    # println("size before", size(vertices))
    normals = normals'
    triangles = triangles'
    vertices = vertices'
    V = V'

    println("triangles: ", size(triangles), " vertices:  ", size(vertices), "V: ", size(V))


    Nvertices = size(vertices, 2)
    gradF = zeros(3,Nvertices)

    for i = 1:5

        mask = findall(x-> x == i, faces)

        num = length(mask)
        println("# kaimu $i-ajam triangle: ", num)
        (row, col) = zeros(Int64, 1, num), zeros(Int64, 1, num)

        for n in 1:num
            row[n], col[n] = mask[n][1], mask[n][2]
        end

        println("rows:", row)
        println("cols", col)
        #(row,col) = find(triangles == i); # atrod virsotnei atbilstoshos trijstuurus

        for this_triang = 1:length(row)
            j1 = (row[this_triang])%(3) + 1 # from 1 to 3
            #println("pre j2: ", (row[this_triang] + 1)%(3))
            j2 = (row[this_triang] + 1)%(3) +1 # from 1 to 3
            #j22 = (row[this_triang] + 1)%(3) +1 # from 1 to 3
            #println("$j2, $j22")
            # radiusvektori uz sho trijsturu virsotneem
            # pirmais atbilst i

            println()
            println("triang $i, kaims = $this_triang ----------------------------------------------------------------------------")

            println()
            println("x index")
            println("trig size = ", size(triangles))
            println("vert size = ", size(vertices))

            println("$(row[this_triang]), $(col[this_triang]) -> ", triangles[row[this_triang], col[this_triang]])
            println("$j1, $(col[this_triang]) -> ", triangles[j1, col[this_triang]])
            println("$j2, $(col[this_triang]) -> ", triangles[j2, col[this_triang]])

            println(size(triangles))
            #println(triangles'[1:20, :])
            println(triangles[1, 1:15])
            println(triangles[2, 1:15])
            println(triangles[3, 1:15])
            println()

            println(size(vertices))
            #
            # x = (vertices[:,triangles[row[this_triang],col[this_triang]]],
            #     vertices[:,triangles[j1, col[this_triang]]],
            #     vertices[:,triangles[j2, col[this_triang]]])

            x = [vertices[:,triangles[row[this_triang],col[this_triang]]],
                vertices[:,triangles[j1, col[this_triang]]],
                vertices[:,triangles[j2, col[this_triang]]]]
            println(x)
            #println()
            v = [V[:, triangles[row[this_triang],col[this_triang]]],
                V[:, triangles[j1, col[this_triang]]],
                V[:, triangles[j2, col[this_triang]]]]

            # this_hsq = [hsq[triangles[row[this_triang],col[this_triang]]]
            #     hsq[triangles[row[this_triang],j1]]
            #     hsq[triangles[row[this_triang],j2]]]

            println("x: ")
            println(typeof(x), size(x))
            println(x[1,:])
            println(x[1])
            println(x[2,:])
            println(x[3,:])

            a = norm(x[2,:] - x[1,:])
            b = norm(x[3,:] - x[2,:])
            c = norm(x[1,:] - x[3,:])

            println("a b c = $a  //  $b  //  $c")

            try
                global Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
                println("first Cdelta = $Cdelta")
            catch
                println()
                println("lalalalalalalal")
                println(" SADNESS    SS")
                global Cdelta = 0
                break
            end

            #(c,b,a) = sort([a,b,c])
            #S = 0.25 * sqrt((a+(b+c)) * (c-(a-b)) * (c+(a-b)) * (a+(b-c)))
            #Cdelta = S / (a^2 + b^2 + c^2)
            #Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

            A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)
            B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)
            C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)

            dCdeltadt = -A * dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) +
                        -B * dot(x[3,:] - x[2,:], v[3,:] - v[2,:]) +
                        -C * dot(x[1,:] - x[3,:], v[1,:] - v[3,:])

            # gradF[i,:] = gradF[i,:] + ...
            #             0.4 / Cdelta^2 * 2 * dCdeltadt*( A*(x[2,:] - x[1,:]) + C*(x[3,:] - x[1,:])) + ...
            #             -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[2])) - 0.5*(this_hsq[1] + this_hsq[2]) / a^4)^2 *(x[2,:] - x[1,:]) + ...
            #             -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[3])) - 0.5*(this_hsq[1] + this_hsq[3]) / c^4)^2 *(x[3,:] - x[1,:])

            #println(sizeof(gradF))
            #println(sizeof(x))
            #println(sizeof(x))
            #println(sizeof(x))
            println("A B C = $A  //  $B  //  $C")
            println("Cdelta = $Cdelta,   dCdt = $dCdeltadt")


            t1 = 0.4 / Cdelta^2 * 2 * dCdeltadt .* ( A.*(x[2,:] - x[1,:]) .+ C.*(x[3,:] - x[1,:]))
            t2 = -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) .* (x[2,:] - x[1,:])
            t3 = -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) .* (x[3,:] - x[1,:])
            println()
            println("$(gradF[:,i]), $(t1[1]), $(t2[1]), $(t3[1])")

            gradF[:,i] = gradF[:,i] + t1[1] + t2[1] + t3[1]



        end
        println(gradF[:,i])
        tang_proj = I - normals[:,i] * normals[:,i]'
        println(tang_proj)
        gradF[:,i] = tang_proj * gradF[:,i]
        println(gradF[:,i])
        println()

    end

    println()
    println(gradF)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    return gradF'

end



function make_gradF(normals,triangles, vertices, V)

    # normals = normals'
    # triangles = triangles'
    # vertices = vertices'
    # V = V'

    Nvertices = size(vertices, 1)

    gradF = zeros(Nvertices,3)

    for i = 1:Nvertices


        mask = findall(x-> x == i, faces)

        num = length(mask)
        (row, col) = zeros(Int64, 1, num), zeros(Int64, 1, num)

        for n in 1:num
            row[n], col[n] = mask[n][1], mask[n][2]
        end
        #(row,col) = find(triangles == i); # atrod virsotnei atbilstoshos trijstuurus

        for this_triang = 1:length(row)
            j1 = (col[this_triang])%(3) + 1; # from 1 to 3
            j2 = (col[this_triang] + 1)%(3) +1; # from 1 to 3

            # radiusvektori uz sho trijsturu virsotneem
            # pirmais atbilst i
            x = [vertices[triangles'[row[this_triang],col[this_triang]],:]
                vertices[triangles'[row[this_triang],j1],:]
                vertices[triangles'[row[this_triang],j2],:]]

            v = [V[triangles'[row[this_triang],col[this_triang]],:]
                V[triangles'[row[this_triang],j1],:]
                V[triangles'[row[this_triang],j2],:]]

            # this_hsq = [hsq[triangles[row[this_triang],col[this_triang]]]
            #     hsq[triangles[row[this_triang],j1]]
            #     hsq[triangles[row[this_triang],j2]]]

            a = norm(x[2,:] - x[1,:])
            b = norm(x[3,:] - x[2,:])
            c = norm(x[1,:] - x[3,:])

            # try
            #     global Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
            # catch
            #     global Cdelta = 0
            # end

            Cdelta = 0.25 * sqrt(1 - 2*(a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)

            A = (a^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)
            B = (b^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)
            C = (c^2 * (a^2 + b^2 + c^2) - a^4 - b^4 - c^4 ) /
                (4*Cdelta * (a^2 + b^2 + c^2)^3)

            dCdeltadt = -A * dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) +
                        -B * dot(x[3,:] - x[2,:], v[3,:] - v[2,:]) +
                        -C * dot(x[1,:] - x[3,:], v[1,:] - v[3,:])

            # gradF[i,:] = gradF[i,:] + ...
            #             0.4 / Cdelta^2 * 2 * dCdeltadt*( A*(x[2,:] - x[1,:]) + C*(x[3,:] - x[1,:])) + ...
            #             -4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[2])) - 0.5*(this_hsq[1] + this_hsq[2]) / a^4)^2 *(x[2,:] - x[1,:]) + ...
            #             -4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (1/(0.5*(this_hsq[1] + this_hsq[3])) - 0.5*(this_hsq[1] + this_hsq[3]) / c^4)^2 *(x[3,:] - x[1,:])
            gradF[i,:] = gradF[i,:] .+
                        0.4 / Cdelta^2 * 2 * dCdeltadt*( A*(x[2,:] - x[1,:]) + C*(x[3,:] - x[1,:])) .+
                        .-4*dot(x[2,:] - x[1,:], v[2,:] - v[1,:]) * (x[2,:] - x[1,:]) .+
                        .-4*dot(x[3,:] - x[1,:], v[3,:] - v[1,:]) * (x[3,:] - x[1,:])


        end


    end

    return gradF

end


function make_Vvecs_conjgrad(normals,triangles, vertices, vvecs, epsilon, maxIters)

# [k1, k2] = principal_curvatures[CDE]; # k1 >= k2
# LAMBDA = k1.^2 + k2.^2 + 0.004
# K = 4/(sqrt(3) * size(triangles,1)) * sum(LAMBDA.^0.25 .* deltaS)
# hsq = K * LAMBDA.^(-0.25)
    triangles = triangles'
    vertices = vertices'
    normals = normals'
    vvecs = vvecs'

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
        println(F)
        println((Fp-F)/F)

        if (Fp - F)/F < epsilon
            break
        end
    end

    return V
end
