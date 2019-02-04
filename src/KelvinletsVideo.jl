module KelvinletsVideo
    include("tetrahedronInterpolator.jl")
    using Images, ProgressMeter, ImageView, LinearAlgebra
    export KelvinletsVideoObject, grab

    struct KelvinletsVideoObject
        a::Float64
        b::Float64
        c::Float64
        sizeX::Int64
        sizeY::Int64
        sizeZ::Int64
        video::AbstractArray{RGB{N0f8}, 3}
        function KelvinletsVideoObject(video::AbstractArray{RGB{N0f8}, 3},
                                        ν::Float64,
                                        μ::Float64
                )::KelvinletsObject

            a = 1 / (4pi * μ)
            b = a / (4(1 - ν))
            c = 2 / (3a- 2b)

            sizeY, sizeX, sizeZ = size(video)

            new(a, b, c, sizeX, sizeY, sizeZ, video)
        end
    end
    function __applyVariation__(object::KelvinletsVideoObject,
                                pressurePoint::Array{Int64, 1},
                                variationCalculator,
                                retardationFunction)::Array{RGB{N0f8}, 3}

        allΔ = zeros(object.sizeY, object.sizeX, object.sizeZ, 3)
        for i=1:object.sizeY
            for j=1:object.sizeX
                for k=1:object.sizeZ
                    Δ = variationCalculator([i, j, k])

                    dx1 = j
                    dx2 = object.sizeX - j
                    dy1 = i
                    dy2 = object.sizeY - i
                    dz1 = k
                    dz2 = object.sizeZ - k

                    dx = min(dx1, dx2)
                    dy = min(dy1, dy2)
                    dz = min(dz1, dz2)

                    y = 2(object.sizeY/2 - dy)/object.sizeY
                    x = 2(object.sizeX/2 - dx)/object.sizeX
                    z = 2(object.sizeZ/2 - dz)/object.sizeZ

                    Δ[1] *= retardationFunction(y)
                    Δ[2] *= retardationFunction(x)
                    Δ[3] *= retardationFunction(z)

                    Δ += [i, j, k]
                    
                    allΔ[i, j, k, :] = Δ[:]

                end
            end
        end
        __interpolateVariation__(object, allΔ)
    end
    function __interpolateVariation__(object::KelvinletsVideoObject,
                                      allΔ::Array{Float64, 4})::Array{RGB{N0f8}, 3}

        interpVideo = fill(RGB{N0f8}(0, 0, 0), object.sizeY, object.sizeX, object.sizeZ)

        rasterize = function(A, B, C, D)
            colorA = object.video[A[1], A[2], A[3]]
            colorB = object.video[B[1], B[2], B[3]]
            colorC = object.video[C[1], C[2], C[3]]
            colorD = object.video[D[1], D[2], D[3]]
            TetrahedronInterpolator.rasterization(interpVideo,
                            allΔ[A[1], A[2], A[3], :],
                            allΔ[B[1], B[2], B[3], :],
                            allΔ[C[1], C[2], C[3], :],
                            allΔ[D[1], D[2], D[3], :],
                            colorA, colorB, colorC, colorD
            )
        end
        @showprogress for i=1:object.sizeY-1
            for j=1:object.sizeX-1
                for k=1:object.sizeZ-1
                    #(0,0,0), (1,0,1), (1,0,0), (1,1,1)
                    rasterize([i, j, k],
                              [i+1, j, k+1],
                              [i+1, j, k],
                              [i+1, j+1, k+1]
                             )
                    #(0,0,0), (1,0,0), (1,1,0), (1,1,1)
                    rasterize([i, j, k],
                              [i+1, j, k],
                              [i+1, j+1, k],
                              [i+1, j+1, k+1]
                             )
                    #(0,0,0), (0,1,0), (1,1,0), (1,1,1)
                    rasterize([i, j, k],
                              [i, j+1, k],
                              [i+1, j+1, k],
                              [i+1, j+1, k+1]
                             )
                    #(0,0,0), (0,1,0), (0,1,1), (1,1,1)
                    rasterize([i, j, k],
                              [i, j+1, k],
                              [i, j+1, k+1],
                              [i+1, j+1, k+1]
                             )
                    #(0,0,0), (0,1,1), (0,0,1), (1,1,1)
                    rasterize([i, j, k],
                              [i, j+1, k+1],
                              [i, j, k+1],
                              [i+1, j+1, k+1]
                             )
                    #(0,0,0), (1,0,1), (0,0,1), (1,1,1)
                    rasterize([i, j, k],
                              [i+1, j, k+1],
                              [i, j, k+1],
                              [i+1, j+1, k+1]
                             )
                end
            end
        end
        interpVideo
    end
    function grab(object::KelvinletsVideoObject,
                  x0::Array{Int64, 1},
                  force::Array{Float64, 1},
                  ϵ::Float64)::Array{RGB{N0f8}, 3}
            grabFunc = function(x)
                r = x - x0
                rLength = norm(r)
                rϵ = sqrt(rLength^2 + ϵ^2)
                kelvinState = (((object.a - object.b)/rϵ) * I +
                                (object.b / rϵ^3) * r * r' +
                                (object.a / 2)*((ϵ^2) / rϵ^3) * I)

                object.c * ϵ * kelvinState * force
            end
            retardationFunc = x -> (cos(π*x) + 1)/2
            __applyVariation__(object, x0, grabFunc, retardationFunc)
        end
end