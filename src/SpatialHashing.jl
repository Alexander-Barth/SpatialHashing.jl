module SpatialHashing

lindex(ind,offsets) = sum((ind .- 1) .* offsets) + 1


function hash(ind,max)
    # https://matthias-research.github.io/pages/publications/tetraederCollision.pdf
    large_primes = (73856093, 19349663, 83492791)[1:length(ind)]
    return abs(reduce(xor,ind .* large_primes)) % (max-1) + 1
end

indices(x,h) = unsafe_trunc.(Int,x ./ h) .+ 1



function spatial_hash!(particles,h,limits,table,num_particles)
    table .= 0
    num_particles .= 0
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1

    # count particles with the same hash
    for i = 1:length(particles)
        ind = indices(particles[i].x,h)
        l = hash(ind,length(table))
        table[l] += 1
    end


    # # partial sums
    # start = 0
    # for i = 1:length(table)
    #     start += table[i];
    #     table[i] = start
    # end
    # table[end] = start

    # partial sums
    for i = 2:length(table)
        table[i] += table[i-1]
    end

    # fill-in

    for i = 1:length(particles)
        ind = indices(particles[i].x,h)
        l = hash(ind,length(table))
        table[l] -= 1
        num_particles[table[l]+1] = i
    end

   # table .+= 1

    return nothing
end

function spatial_hash(particles,h,limits)
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1
    table = zeros(Int,prod(sz)+1)
    num_particles = zeros(Int,length(particles))
    @time spatial_hash!(particles,h,limits,table,num_particles)
    return (; table, num_particles, h, sz)
end

function each_near(fun,x,search_range,spatial_index,visited)
    visited .= 0
    table,num_particles,h,sz = spatial_index
    N = length(x)
    ind = indices(x,h)

    search = ntuple(N) do i
        max(1,(ind[i]-search_range)):min(sz[i],(ind[i]+search_range))
    end

    #@show search,sz
    for ind2 in CartesianIndices(search)
        l = hash(Tuple(ind2),length(table))
        #@show ind2,sz,l,length(table)
        for i = table[l]:(table[l+1]-1)
            # must inline to avoid allocation

            j = num_particles[i+1]

            if visited[j] == 0
                @inline fun(j)
                visited[j] = 1
            end
        end
    end
end

end
