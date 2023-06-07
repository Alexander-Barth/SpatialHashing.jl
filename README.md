# SpatialHashing

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Alexander-Barth.github.io/SpatialHashing.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Alexander-Barth.github.io/SpatialHashing.jl/dev/)
[![Build Status](https://github.com/Alexander-Barth/SpatialHashing.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Alexander-Barth/SpatialHashing.jl/actions/workflows/CI.yml?query=branch%3Amain)


From a cloud of points in an n-dimensional space,
this package allows to identify all points near a given search points.
The search is done using spatial hashing requiring `O(log(N))` operations per search
where `N` is the number of points.

The code is based on:
* Matthias Müller: [Blazing Fast Neighbor Search with Spatial Hashing](https://matthias-research.github.io/pages/tenMinutePhysics/11-hashing.pdf), Ten Minute Physics, 2022
* Matthias Teschner, Bruno Heidelberger, Matthias Müller, Danat Pomerantes, Markus H. Gross: [Optimized Spatial Hashing for Collision Detection of Deformable Objects](https://matthias-research.github.io/pages/publications/tetraederCollision.pdf). VMV 2003: 47-54


The implementation is intended to use minimal dependencies (currently there are no dependencies) and compilable using the julia subset allowed by [GPUCompiler.jl](https://github.com/JuliaGPU/GPUCompiler.jl).

Documentation is available [here](https://alexander-barth.github.io/SpatialHashing.jl/dev/).
