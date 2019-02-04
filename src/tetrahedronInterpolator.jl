module TetrahedronInterpolator
    export rasterization
    using Images
    function planeEquation(a::Array{Float64, 1}, b::Array{Float64, 1}, c::Array{Float64, 1})::Array{Float64, 1}
         @fastmath @inbounds ab::Array{Float64, 1} = [b[1] - a[1], b[2] - a[2], b[3] - a[3]]
         @fastmath @inbounds ac::Array{Float64, 1} = [c[1] - a[1], c[2] - a[2], c[3] - a[3]]
        normal::Array{Float64, 1} = cross(ab, ac)
         @fastmath @inbounds d::Float64 = - (normal[1]*a[1] + normal[2]*a[2] + normal[3]*a[3])
         [normal[1], normal[2], normal[3], d]
    end

    function getEquations(A::Array{Float64, 1}, B::Array{Float64, 1},
                          C::Array{Float64, 1}, D::Array{Float64, 1})::Array{Array{Float64, 1}, 1}
        eq1::Array{Float64, 1} = planeEquation(A, B, C)
        if (eq1[1:3]' * D + eq1[4])[1] < 0
            eq1 = planeEquation(A, C, B)
        end
        eq2::Array{Float64, 1} = planeEquation(A, B, D)
        if (eq2[1:3]' * C + eq2[4])[1] < 0
            eq2 = planeEquation(A, D, B)
        end
        eq3::Array{Float64, 1} = planeEquation(A, C, D)
        if (eq3[1:3]' * B + eq3[4])[1] < 0
            eq3 = planeEquation(A, D, C)
        end
        eq4::Array{Float64, 1} = planeEquation(B, C, D)
        if (eq4[1:3]' * A + eq4[4])[1] < 0
            eq4 = planeEquation(B, D, C)
        end
        [eq1, eq2, eq3, eq4]
    end

    function checkPoint(eq1::Array{Float64, 1}, eq2::Array{Float64, 1},
                        eq3::Array{Float64, 1}, eq4::Array{Float64, 1},
                        P::Array{Int64, 1})::Bool
        (eq1[1:3]' * P + eq1[4]) >= 0 &&
        (eq2[1:3]' * P + eq2[4]) >= 0 &&
        (eq3[1:3]' * P + eq3[4]) >= 0 &&
        (eq4[1:3]' * P + eq4[4]) >= 0
    end

    function pointToPlane(point::Array{Float64, 1}, a::Array{Float64, 1},
                          b::Array{Float64, 1}, c::Array{Float64, 1})::Float64
        abc::Array{Float64, 1} = planeEquation(a, b, c)

        b::Float64 = abc[1]
        a::Float64 = abc[2]
        c::Float64 = abc[3]
        d::Float64 = abc[4]

         α::Float64 = - (a*point[2] + b*point[1] + c*point[3] + d)/(a^2 + b^2 + c^2)

        pointOnPlane::Array{Float64, 1} = α*abc[1:3] + point

        norm(point - pointOnPlane)
    end

    function interpolateColors(P::Array{Int64, 1},
                               A::Array{Float64, 1}, B::Array{Float64, 1},
                               C::Array{Float64, 1}, D::Array{Float64, 1},
                               colorA::RGB{N0f8}, colorB::RGB{N0f8},
                               colorC::RGB{N0f8}, colorD::RGB{N0f8})::RGB{N0f8}
        influenceA::Float64 = pointToPlane(Array{Float64}(P), B, C, D) / pointToPlane(A, B, C, D)
        influenceB::Float64 = pointToPlane(Array{Float64}(P), A, C, D) / pointToPlane(B, A, C, D)
        influenceC::Float64 = pointToPlane(Array{Float64}(P), A, B, D) / pointToPlane(C, A, B, D)
        influenceD::Float64 = pointToPlane(Array{Float64}(P), A, B, C) / pointToPlane(D, A, B, C)

        influenceA * colorA +
        influenceB * colorB +
        influenceC * colorC +
        influenceD * colorD
    end

    function setBBOX(A::Array{Float64, 1}, B::Array{Float64, 1}, C::Array{Float64, 1}, D::Array{Float64, 1},
                    Δx::Int64, Δy::Int64, Δz::Int64)::Array{Real, 1}
        maxY::Real = max(A[1], B[1], C[1], D[1])
        minY::Real = min(A[1], B[1], C[1], D[1])

        if maxY > Δy
            maxY = Δy
        end
        if minY < 1
            minY = 1
        end

        maxX::Real = max(A[2], B[2], C[2], D[2])
        minX::Real = min(A[2], B[2], C[2], D[2])
        if maxX > Δx
            maxX = Δx
        end
        if minX < 1
            minX = 1
        end

        maxZ::Real = max(A[3], B[3], C[3], D[3])
        minZ::Real = min(A[3], B[3], C[3], D[3])
        if maxZ > Δz
            maxZ = Δz
        end
        if minZ < 1
            minZ = 1
        end
        [maxY, minY, maxX, minX, maxZ, minZ]
    end

    function checkAndSet(eq1::Array{Float64, 1}, eq2::Array{Float64, 1},
                         eq3::Array{Float64, 1}, eq4::Array{Float64, 1},
                         video::Array{RGB{N0f8},3}, A::Array{Float64, 1}, B::Array{Float64, 1},
                         C::Array{Float64, 1}, D::Array{Float64, 1}, P::Array{Int64, 1},
                         colorA::RGB{N0f8}, colorB::RGB{N0f8}, colorC::RGB{N0f8}, colorD::RGB{N0f8})
        if checkPoint(eq1, eq2, eq3, eq4, P)
            @inbounds video[P[1], P[2], P[3]] = interpolateColors(P, A, B, C, D, colorA, colorB, colorC, colorD)
        end
    end

    function rasterization(video::Array{RGB{N0f8},3},
                           A::Array{Float64, 1}, B::Array{Float64, 1},
                           C::Array{Float64, 1}, D::Array{Float64, 1},
                           colorA::RGB{N0f8}, colorB::RGB{N0f8},
                           colorC::RGB{N0f8}, colorD::RGB{N0f8})::Array{RGB{N0f8}, 3}
        maxY::Float64, minY::Float64, maxX::Float64, minX::Float64, maxZ::Float64, minZ::Float64 = setBBOX(A, B, C, D,
                                                                    size(video)[2], size(video)[1], size(video)[3])
        eq1, eq2, eq3, eq4 = getEquations(A, B, C, D)
        for i=floor(Int64, minY):ceil(Int64, maxY)
            for j=floor(Int64, minX):ceil(Int64, maxX)
                @simd for k=floor(Int64, minZ):ceil(Int64, maxZ)
                     checkAndSet(eq1, eq2, eq3, eq4, video, A, B, C, D, [i, j, k], colorA, colorB, colorC, colorD)
                end
            end
        end
        video
    end
end
