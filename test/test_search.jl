using SpatialHashing: spatial_hash!, spatial_hash, each_near
using SpatialHashing
using Test
using StaticArrays
using LinearAlgebra
using StableRNGs


N = 2

for nparticles in (10,10000)
   for h = (0.01, 0.4,0.8)

#h = 0.8
#nparticles = 10000
        @show h,nparticles

        rng = StableRNG(42)
        nparticles = 10000
        xpos = rand(N,nparticles)
        particles = [(x=SVector(tuple(rand(rng,N)...)),) for i = 1:nparticles]

        limits = (2,2)


        visited = zeros(Int,nparticles)
        sz = unsafe_trunc.(Int,limits ./ h) .+ 1
        table = zeros(Int,prod(sz)+1)
        num_particles = zeros(Int,length(particles))
        limits = Tuple(limits)
        @time spatial_hash!(particles,h,limits,table,num_particles)


        #@code_warntype spatial_hash!(particles,h,limits,table,num_particles)



        # #locations(i) = particles[i].x
        # using StaticArrays
        # import Base: size, getindex

        # struct Location3{N,T,TC <: AbstractVector{Particle{N,T}}} <: AbstractVector{SVector{N,T}}
        #     p::TC
        # end
        # @inline Base.size(l::Location3) = size(l.p)
        # @inline Base.@propagate_inbounds Base.getindex(l::Location3,index) = l.p[index].x


        # l = Location3(particles);

        spatial_index = spatial_hash(particles,h,limits)

        #@test table[end] == length(particles)


        search_range = 1
        radius2 = h^2
        offset = first(CartesianIndices(ntuple(i -> -search_range:search_range,N)))


        function find_near!(spatial_index,particles,x,search_range,r2max,near_indices,visited)
            nfound = 0

            @inline each_near(x,search_range,spatial_index,visited) do j
                pj = particles[j]
	            rij = pj.x - x
	            r2 = norm(rij)^2

	            if r2 < r2max
                    nfound += 1
                    #            if nfound <= length(near_indices)
                    #nfound += 1
                    @inbounds near_indices[nfound] = j
                    #            end
                end
            end
            return nfound
        end

        near_indices = zeros(Int,length(particles))
        r2max = h^2;

        for i = 1:length(particles)
            #@show i
            pi = particles[i]
            x = pi.x

            nfound = find_near!(spatial_index,particles,x,search_range,r2max,near_indices,visited)

            near = near_indices[1:nfound]

            if i < 100
                near_ref = Int[]

                for j = 1:length(particles)
                    pj = particles[j]
	                rij = pj.x - pi.x
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
