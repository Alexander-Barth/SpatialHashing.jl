module SpatialHashing

lindex(ind,offsets) = sum((ind .- 1) .* offsets) + 1


function hash(ind,max)
    # The 3 prime numbers are from:
    # https://matthias-research.github.io/pages/publications/tetraederCollision.pdf
    large_primes = (73856093, 19349663, 83492791, 22335757, 98746231, 10000019)[1:length(ind)]
    return abs(reduce(xor,ind .* large_primes)) % (max-1) + 1
end

indices(x,h) = unsafe_trunc.(Int,x ./ h) .+ 1

"""
    SpatialHashing.spatial_hash!(particles,h,limits,table,num_particles)


Initialize the data structure for spatial hashing.
"""
function spatial_hash!(particles,h,limits,table,num_particles)
    table .= 0
    num_particles .= 0
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1

    # count particles with the same hash
    for i = 1:length(particles)
        ind = indices(particles[i],h)
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
        ind = indices(particles[i],h)
        l = hash(ind,length(table))
        table[l] -= 1
        num_particles[table[l]+1] = i
    end

   # table .+= 1

    return nothing
end


"""
    index = SpatialHashing.spatial_hash(particles,h,limits)


Initialize the data structure `index` for spatial hashing using the vector
of position `particles` and resolution `h`. The i-th coordinate is bounded by
`0` and `limits[i]`. The list `particles` is a iterable of coordinates such that
the position of the j-th particle is given by `particles[j]`.

The `index` allows to localize all points near a given reference point in
`O(log(n))` operations where `n` is the number of particles.

The hash function is currently setup for 1 to 6-dimensional spaces.

The parameter `h` is used to define a grid of cell. Each particle belongs to
one grid cell. For every cell a
hash number is computed and the `index` keeps track of all `particles`
with the same hash.

```
   <--h-->
   +------+------+------+------+------+------+  (limits[1],limits[2])
   | ①    |      |      |      |      |      |
   |      |      |      |      |      |      |
   +------+------+------+------+------+------+
   |      |      |      |      |   ⑤  |      |
   |      |      |      |      | ③    |      |
   +------+------+------+------+------+------+
   |      |      |      |      |      |      |
   |      |  ②   |      |      |      |      |
   +------+------+------+------+------+------+
   |      |      |      |  ④   |      |      |
   |      |      |      |      |      |      |
   +------+------+------+------+------+------+
 (0,0)
```
"""
function spatial_hash(particles,h,limits)
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1
    table = zeros(Int,prod(sz)+1)
    num_particles = zeros(Int,length(particles))
    spatial_hash!(particles,h,limits,table,num_particles)
    return (; table, num_particles, h, sz)
end

"""
    SpatialHashing.each_near(fun,x,search_range,spatial_index,visited)

Evaluates `fun` for every index `j` in `particles` which is near `x`.
The search is limited to ±spatial_index grid cells.
The vector `visited` is an vector of booleans of the size `particles`
to keep track which particles have been found so far.

It is important that the function `fun` checks that the proposed point with
index `j` is really close.

The following is a complete example in 2 dimensions:

```julia
nparticles = 10000
particles = [Tuple(rand(2)) for i = 1:nparticles]
h = 0.1
limits = (1,1)
index = SpatialHashing.spatial_hash(particles,h,limits)
x = (0.2,0.2)
r2max = 0.1^2
search_range = 1
visited = falses(nparticles)
SpatialHashing.each_near(x,search_range,spatial_index,visited) do j
   r2 = sum((particles[j] .- x).^2)
   if r2 < r2max
      println("point with index \$j is near (distance is \$(sqrt(r2)))")
   end
end
```

"""
function each_near(fun,x,search_range,spatial_index,visited)
    visited .= 0
    table,num_particles,h,sz = spatial_index
    N = length(x)
    ind = indices(x,h)

    search = ntuple(N) do i
        max(1,(ind[i]-search_range)):min(sz[i],(ind[i]+search_range))
    end

    for ind2 in CartesianIndices(search)
        l = hash(Tuple(ind2),length(table))
        for i = table[l]:(table[l+1]-1)
            j = num_particles[i+1]

            if visited[j] == 0
                # must inline to avoid allocation
                @inline fun(j)
                visited[j] = 1
            end
        end
    end
end

end
