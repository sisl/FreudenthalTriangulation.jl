# Usage

## Getting Started

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
