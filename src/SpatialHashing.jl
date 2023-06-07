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
    SpatialHashing.spatial_hash!(points,h,limits,table,num_points)


Initialize the data structure for spatial hashing in-place.
`table` is a vector of integers with the length of
`prod(trunc.(Int,limits ./ h) .+ 1)+1`, `num_points` is a vector
of integers with the same length as `points`.

See `SpatialHashing.spatial_hash` for explanations of the other parameters.
"""
function spatial_hash!(points,h,limits,table,num_points)
    table .= 0
    num_points .= 0
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1

    # count points with the same hash
    for i = 1:length(points)
        ind = indices(points[i],h)
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

    for i = 1:length(points)
        ind = indices(points[i],h)
        l = hash(ind,length(table))
        table[l] -= 1
        num_points[table[l]+1] = i
    end

   # table .+= 1

    return nothing
end


"""
    spatial_index = SpatialHashing.spatial_hash(points,h,limits)


Initialize the data structure `spatial_index` for spatial hashing using the vector
of position `points` and resolution `h`. The i-th coordinate is bounded by
`0` and `limits[i]`. The list `points` is a iterable of coordinates such that
the position of the j-th particle is given by `points[j]`.

The `index` allows to localize all points near a given reference point in
`O(log(n))` operations where `n` is the number of points.

The hash function is currently setup for 1- up to 6-dimensional spaces.

The parameter `h` is used to define a grid of cell. Each particle belongs to
one grid cell. For every cell a
hash number is computed and the `index` keeps track of all `points`
with the same hash.

```
   <──h──>                           (limits[1],limits[2])
   ┏━━━━━━┯━━━━━━┯━━━━━━┯━━━━━━┯━━━━━━┯━━━━━━┓
   ┃      │      │      │      │      │      ┃
   ┃      │      │      │      │      │      ┃
   ┠──────┼──────┼──────┼──────┼──────┼──────┨
   ┃  ①   │      │      │      │      │      ┃
   ┃      │      │      │      │      │      ┃
   ┠──────┼──────┼──────┼──────┼──────┼──────┨
   ┃      │      │      │  ④   │      │      ┃
   ┃      │      │      │      │      │      ┃
   ┠──────┼──────┼──────┼──────┼──────┼──────┨
   ┃      │  ②   │      │      │      │      ┃
   ┃      │      │      │      │      │      ┃
   ┠──────┼──────┼──────┼──────┼──────┼──────┨
   ┃      │      │      │      │    ⑤ │      ┃
   ┃      │      │      │      │ ③    │      ┃
   ┗━━━━━━┷━━━━━━┷━━━━━━┷━━━━━━┷━━━━━━┷━━━━━━┛
 (0,0)
```
"""
function spatial_hash(points,h,limits)
    sz = unsafe_trunc.(Int,limits ./ h) .+ 1
    table = zeros(Int,prod(sz)+1)
    num_points = zeros(Int,length(points))
    spatial_hash!(points,h,limits,table,num_points)
    return (; table, num_points, h, sz, limits)
end

"""
    SpatialHashing.update!(spatial_index,points)

Update the spatial `spatial_index` using the new positions of
the `points`. The number of points is assumed to be the same
as during initialization with `SpatialHashing.spatial_hash`.
"""
function update!(spatial_index,points)
    table,num_points,h,sz,limits = spatial_index
    spatial_hash!(points,h,limits,table,num_points)
    return (; table, num_points, h, sz, limits)
end

"""
    SpatialHashing.each_near(fun,x,search_range,spatial_index,visited)

Evaluates `fun` for every index `j` in `points` which is near `x`.
The search is limited to ±spatial_index grid cells.
The vector `visited` is an vector of booleans of the size `points`
to keep track which points have been found so far.

It is important that the function `fun` checks that the proposed point with
index `j` is actually close as different grid cells may have the same hash
(hash-collision).

The following is a complete example in 2 dimensions:

```julia
npoints = 10000
points = [Tuple(rand(2)) for i = 1:npoints]
h = 0.1
limits = (1,1)
index = SpatialHashing.spatial_hash(points,h,limits)
x = (0.2,0.2)
r2max = 0.1^2
search_range = 1
visited = falses(npoints)
SpatialHashing.each_near(x,search_range,spatial_index,visited) do j
   r2 = sum((points[j] .- x).^2)
   if r2 < r2max
      println("point with index \$j is near (distance is \$(sqrt(r2)))")
   end
end
```

"""
function each_near(fun,x,search_range,spatial_index,visited)
    visited .= 0
    table,num_points,h,sz,limits = spatial_index
    N = length(x)
    ind = indices(x,h)

    search = ntuple(N) do i
        max(1,(ind[i]-search_range)):min(sz[i],(ind[i]+search_range))
    end

    for ind2 in CartesianIndices(search)
        l = hash(Tuple(ind2),length(table))
        for i = table[l]:(table[l+1]-1)
            j = num_points[i+1]

            if visited[j] == 0
                # must inline to avoid allocation
                @inline fun(j)
                visited[j] = 1
            end
        end
    end
end

end
