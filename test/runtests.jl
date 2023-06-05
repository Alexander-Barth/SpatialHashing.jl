using SpatialHashing: spatial_hash!, spatial_hash, each_near
using SpatialHashing
using Test
using StaticArrays
using LinearAlgebra
using StableRNGs

@testset "SpatialHashing.jl" begin
    include("test_search.jl")
end
