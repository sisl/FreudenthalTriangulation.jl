# FreudenthalTriangulations.jl

[![Build Status](https://travis-ci.org/sisl/FreudenthalTriangulations.jl.svg?branch=main)](https://travis-ci.org/sisl/FreudenthalTriangulations.jl)
[![Coverage](https://coveralls.io/repos/github/sisl/FreudenthalTriangulations.jl/badge.svg?branch=main)](https://coveralls.io/github/sisl/FreudenthalTriangulations.jl?branch=main)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sisl.github.io/FreudenthalTriangulations.jl/)

This package allows users to use Freudenthal triangulation for [functional value approximation](https://sisl.github.io/FreudenthalTriangulation.jl/concepts/#Concepts).

Specifically it allows the user to find the simplex of a point `x` in Freudenthal space and the barycentric coordinates of that point with respect to the its simplex.

## Installation

Start Julia and run the following command:

```julia
Pkg.add("FreudenthalTriangulations")
```

## Usage

To use the FreudenthalTriangulation module, begin your code with

```julia
using FreudenthalTriangulations
```
<!--
## Finding Simplex and Barycentric Coordinates

To find the vertices of the simplex around a point `x` in Freudenthal space, run
```julia
V = freudenthal_simplex(x)
```
Then to find barycentric coordinates of `x` with respect to the simplex, run
```julia
coords = barycentric_coordinates(x, V)
```
Note that these vertices must be in the same order as provided by `freudenthal_simplex`.

To calculate the simplex and the barycentric coordinates with one function, run
```julia
V, coords = freudenthal_simplex_and_coords(x)
```
where `V` is the simplex vertices in the Freudenthal space and `coords` will be the barycentric coordinates with respect to the simplex.

To calculate the simplex and barycentric coordinates without allocating more memory, run
```julia
V, coords = freudenthal_simplex_and_coords!(x, V, coords)
```
so that `V` will be filled with the simplex vertices in the Freudenthal space and `coords` will be filled with the associated coordinates.

For these functions the requirements are
- `x::Vector{Float64}` The point in Freudenthal space
- `V::Vector{Vector{Float64}}` The vertices of the simplex around `x` in Freudenthal space
- `coords::Vector{Float64}` The barycentric coordinates of `x` with respect to the simplex -->

## Credits

Contributors to this package include Sidhart Krishnan, Tim Wheeler, and Mykel Kochenderfer.
