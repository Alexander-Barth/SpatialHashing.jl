module SpatialHashing
import Base: iterate

function hash(ind,max)
    # The 3 first prime numbers are from:
    # https://matthias-research.github.io/pages/publications/tetraederCollision.pdf
    large_primes = (73856093, 19349663, 83492791, 22335757, 98746231, 10000019)[1:length(ind)]
    return abs(reduce(xor,ind .* large_primes)) % (max-1) + 1
end

indices(x,h) = unsafe_trunc.(Int,x ./ h) .+ 1

"""
    SpatialHashing.spatial_hash!(points,h,table,num_points)


Initialize the data structure for spatial hashing in-place.
`table` is a vector of integers with hashes, `num_points` is a vector
of integers with the same length as `points`.

See `SpatialHashing.spatial_hash` for explanations of the other parameters.
"""
function spatial_hash!(points,h,table,num_points)
    table .= 0
    num_points .= 0

    # count points with the same hash
    for i = 1:length(points)
        ind = indices(points[i],h)
        l = hash(ind,length(table))
        table[l] += 1
    end

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

    return nothing
end


"""
    spatial_index = SpatialHashing.spatial_hash(points,h,max_table)


Initialize the data structure `spatial_index` for spatial hashing using the vector
of position `points` and grid with resolution `h`. `max_table` is the length
of the hash table (while size can be tuned to optimize the performance).

The `index` allows to localize all points near a given reference point in
`O(log(n))` operations where `n` is the number of points.

The hash function is currently setup for 1- up to 6-dimensional spaces.

The parameter `h` is used to define a grid of cell. Each particle belongs to
one grid cell. For every cell a
hash number is computed and the `index` keeps track of all `points`
with the same hash.

```
   <──h──>
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
function spatial_hash(points,h,max_table)
    table = zeros(Int,max_table)
    num_points = zeros(Int,length(points))
    spatial_hash!(points,h,table,num_points)
    return (; table, num_points, h)
end

"""
    SpatialHashing.update!(spatial_index,points)

Update the spatial `spatial_index` using the new positions of
the `points`. The number of points is assumed to be the same
as during initialization with `SpatialHashing.spatial_hash`.
"""
function update!(spatial_index,points)
    table,num_points,h = spatial_index
    spatial_hash!(points,h,table,num_points)
    return (; table, num_points, h)
end

"""
    SpatialHashing.each_near(fun,x,search_range,spatial_index,visited)

Evaluates `fun` for every index `j` in `points` which is near `x`.
The search is limited to ± `spatial_index` around the grid cell containing `x`.
The vector `visited` is an vector of booleans of the size `points`
to keep track which points have been found so far.

It is important that the function `fun` checks that the proposed point with
index `j` is actually close as different grid cells may have the same hash
(hash-collision).

The following is a complete example in 2 dimensions:

```julia
npoints = 200
points = [Tuple(rand(2)) for i = 1:npoints]
h = 0.1
max_table = 50
spatial_index = SpatialHashing.spatial_hash(points,h,max_table)
x = (0.2,0.2)
r2max = 0.1^2
search_range = 1
visited = falses(npoints)
SpatialHashing.each_near(x,search_range,spatial_index,visited) do j
   r2 = sum((points[j] .- x).^2)
   if r2 < r2max
      println("point with index ",j," is near (distance is ",sqrt(r2),")")
   end
end

# the same indices as naive search

filter(i -> sum((points[i] .- x).^2) < r2max,1:length(points))
```

"""
function each_near(fun,x,search_range,spatial_index,visited)
    visited .= 0
    table,num_points,h = spatial_index
    N = length(x)
    ind = indices(x,h)

    search = ntuple(N) do i
        (ind[i]-search_range):(ind[i]+search_range)
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

struct Iter{Tindex,Tx,T,Tv,Tc <: CartesianIndices}
    spatial_index::Tindex
    x::Tx
    search_range::T
    visited::Tv
    indices::Tc
end

function init_(it)
    state1 = nothing
    inner = false
    next2 = nothing
    iter2 = nothing
    iii = 0
    l = 0
    next1 = iterate(it.indices)

    if next1 == nothing
        return nothing
    else
        item1,state1 = next1
        (item1,(state1,inner,next2,iter2,iii,l,next1))
    end
end

function next_(it,state)
    (state1,inner,next2,iter2,iii,l,next1) = state
    if inner
        iii += 1
    else
        next1 = iterate(it.indices, state1)
    end


    if next1 == nothing
        return nothing
    else
        item1,state1 = next1
        return (item1,(state1,inner,next2,iter2,iii,l,next1))
    end
end

function iterate(it,state__ = nothing)
    if state__ == nothing
        state1=inner=next2=iter2=iii=l=next1 = nothing
        next = init_(it)
    else
        (next,state1,inner,next2,iter2,iii,l,next1) = state__
        state_ = (state1,inner,next2,iter2,iii,l,next1)
        next = next_(it,state_)
    end

    visited = it.visited

    while next !== nothing
        table,num_points,h = it.spatial_index
        ind2,state_ = next
        (state1,inner,next2,iter2,iii,l,next1) = state_

        if !inner
            (ind2, state1) = next1
            l = hash(Tuple(ind2),length(table))
            iii = table[l]
        end

        if iii > table[l+1]-1
            inner = false
        else
            inner = true
            j = num_points[iii+1]

            if visited[j] == 0
                visited[j] = 1
                return (j,(next,state1,inner,next2,iter2,iii,l,next1))
            end
        end

        state_ = (state1,inner,next2,iter2,iii,l,next1)
        next = next_(it,state_)
    end

    return nothing
end


function inrange!(spatial_index,x,search_range,visited)
    visited .= 0
    table,num_points,h = spatial_index
    N = length(x)
    ind = indices(x,h)

    search = ntuple(N) do i
        (ind[i]-search_range):(ind[i]+search_range)
    end

    return Iter(spatial_index,x,search_range,visited,CartesianIndices(search))
end


end
