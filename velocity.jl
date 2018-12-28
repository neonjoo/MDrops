function InterfaceSpeedZinchenko(points,faces,forcen,etaP,gammap)

    ### A possible improvement would be the interpolation of curvature and normal
    normals = Array{Float64}(undef,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
    
    #phi = Array(Float64,size(points,2))
    velocityn = zeros(Float64,size(points,2))
    
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]
        fx = forcen[xkey]
                
        s = 0
        for ykey in 1:size(points,2)
            if ykey==xkey
                continue
            end

            y = points[:,ykey]
            ny = normals[:,ykey]
            fy = forcen[ykey]

            ### I will need to check a missing 2
            s += vareas[ykey]*1 ./8/pi/etaP* dot(y-x,nx+ny)/norm(y-x)^3*(1-3*dot(y-x,nx)*dot(y-x,ny)/norm(y-x)^2) * gammap

            ### ????????
            s += vareas[ykey]*1 ./8/pi/etaP* ( dot(nx,ny)/norm(x-y) + dot(nx,x -y)*dot(ny,x-y)/norm(x-y)^3 )*(fy - fx)
        end

        velocityn[xkey] = s
    end

    return velocityn
end

