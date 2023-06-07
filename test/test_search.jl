using SpatialHashing: spatial_hash!, spatial_hash, each_near
using SpatialHashing
using Test
using LinearAlgebra
using StableRNGs


function find_near!(spatial_index,points,x,search_range,r2max,near_indices,visited)
    nfound = 0

    @inline each_near(x,search_range,spatial_index,visited) do j
        pj = points[j]
	    rij = pj .- x
	    r2 = norm(rij)^2

	    if r2 < r2max
            nfound += 1
            @inbounds near_indices[nfound] = j
        end
    end
    return nfound
end


N = 2

for npoints in (10,1000)
    for h = (0.01, 0.4, 0.8)
        rng = StableRNG(42)
        points = [tuple(rand(rng,N)...) for i = 1:npoints]

        limits = (2,2)

        visited = zeros(Int,npoints)

        spatial_index = spatial_hash(points,h,limits)

        search_range = 1
        near_indices = zeros(Int,length(points))
        r2max = h^2;

        for i = 1:length(points)
            x = points[i]

            nfound = find_near!(spatial_index,points,x,search_range,r2max,near_indices,visited)

            near = near_indices[1:nfound]

            if i < 100
                near_ref = Int[]

                for j = 1:length(points)
                    pj = points[j]
	                rij = pj .- x
	                r2 = norm(rij)^2

	                if r2 < r2max
                        push!(near_ref,j)
                    end
                end

                @test Set(near_ref) == Set(near)
            end
        end
    end
end
