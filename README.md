# FreudenthalTriangulation

[![Build Status](https://travis-ci.com/SidhartK/FreudenthalTriangulation.jl.svg?branch=master)](https://travis-ci.com/SidhartK/FreudenthalTriangulation.jl)
[![Coverage](https://codecov.io/gh/SidhartK/FreudenthalTriangulation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SidhartK/FreudenthalTriangulation.jl)
[![Coverage](https://coveralls.io/repos/github/SidhartK/FreudenthalTriangulation.jl/badge.svg?branch=master)](https://coveralls.io/github/SidhartK/FreudenthalTriangulation.jl?branch=master)

This package allows users to use Freudenthal triangulation for local value function approximation over a set of discrete belief points.

Specifically it allows the user to find the simplex of a point `x` in Freudenthal space and the barycentric coordinates of that point with respect to the its simplex. It also allows conversion between belief space and Freudenthal space.

## Installation

Start Julia and run the following command:

```julia
Pkg.add("FreudenthalTriangulation")
```

## Usage

To use the FreudenthalTriangulation module, begin your code with

```julia
using FreudenthalTriangulation
```

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
V, coords = freudenthal_simplex(x, V, coords)
```
so that `V` will be filled with the simplex vertices in the Freudenthal space and `coords` will be filled with the associated coordinates.

For these functions the requirements are
- `x::Vector{Float64}` The point in Freudenthal space
- `V::Vector{Vector{Float64}}` The vertices of the simplex around `x` in Freudenthal space
- `coords::Vector{Float64}` The barycentric coordinates of `x` with respect to the simplex

## Moving between Belief and Freudenthal Space

To go from a point `x` in Freudenthal space to a point in belief space use the function
```julia
b = to_belief(x, m)
```
where `m` is the granularity of the triangulation.

To go from a point in belief space `b` to a point in Freudenthal space use the function
```julia
x = to_freudenthal(b, m)
```
where  `m` is again the granularity.

To convert a batch of belief points `B` in a Freudenthal space with granularity `m`, use
```julia
X = to_freudenthal_batch(B, m)
```
For these functions the requirements are
- `m::Int64` Granularity of the Freudenthal triangulation
- `B::AbstractArray` The columns are the belief state points

## Credits

Contributors to this package include Mykel Kochenderfer and Sidhart Krishnan
