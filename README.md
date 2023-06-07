# SpatialHashing

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Alexander-Barth.github.io/SpatialHashing.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Alexander-Barth.github.io/SpatialHashing.jl/dev/)
[![Build Status](https://github.com/Alexander-Barth/SpatialHashing.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Alexander-Barth/SpatialHashing.jl/actions/workflows/CI.yml?query=branch%3Amain)


From a cloud of point in n-dimensional space,
this module allows to identify all points near the given search points.
The search is done using a spatial hashing requiring `O(log(N))` operations per search
where `N` is the number of points.

The code is based on:
* Matthias Müller: [Blazing Fast Neighbor Search with Spatial Hashing](https://matthias-research.github.io/pages/tenMinutePhysics/11-hashing.pdf), Ten Minute Physics, 2022
* Matthias Teschner, Bruno Heidelberger, Matthias Müller, Danat Pomerantes, Markus H. Gross: [Optimized Spatial Hashing for Collision Detection of Deformable Objects](https://matthias-research.github.io/pages/publications/tetraederCollision.pdf). VMV 2003: 47-54
