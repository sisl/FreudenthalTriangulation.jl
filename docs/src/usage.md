# Usage

## Getting Started

To use the FreudenthalTriangulations module, begin your code with

```julia
using FreudenthalTriangulations
```

## FreudenthalTriangulation
The `FreudenthalTriangulation` struct can be used to specify the dimension of the space which is represented by `n` and the `m` represents the granularity of the triangulation.

## Finding Freudenthal Vertices
To find a Freudenthal triangulation of a space of dimension `n` with granularity `m` use the `vertices` function which takes in a `FreudenthalTriangulation` struct and returns a list of vertices. Each of these vertices are represented by `n` dimensional vectors. Thus if we set `T = FreudenthalTriangulation(n, m)` and
```julia
V = vertices(T)
```
then `V` is a list of vertices comprising a freudenthal triangulation of an `n` dimensional space.

## Finding Simplex and Barycentric Coordinates

To find the vertices of the simplex around a point `x` in Freudenthal space and the barycentric coordinates, run
```julia
V, coords = simplex(x)
```
Then we have that `V` is the simplex vertices in the Freudenthal space and `coords` will be the barycentric coordinates with respect to the simplex.

For these functions, the requirements are
- `x::Vector{Float64}` The point in Freudenthal space
- `V::Vector{Vector{Float64}}` The vertices of the simplex around `x` in Freudenthal space
- `coords::Vector{Float64}` The barycentric coordinates of `x` with respect to the simplex

## Operating in Belief Space
This package provides two functions to operate in belief space. First to convert a freudenthal triangulation into belief space use the function
```julia
bv = belief_vertices(T)
```
We see that `bv` is a set of vertices in belief space corresponding to a freudenthal triangulation in Freudenthal space. Given any point in belief space we can find a simplex of points in these belief vertices `bv` and from there approximate the value function at this belief.

This leads us into the second function over the belief space, `belief_simplex`. This function allows us to calculate the belief simplex of a belief.
```julia
B, coords = belief_simplex(T, b)
```
where `T` is a FreudenthalTriangulation struct and `b` is a belief. Then we have that `B` is a vector containing a set of beliefs corresponding to the belief simplex and `coords` is the barycentric coordinates of the belief in the belief simplex.

For these functions, the requirements are
- `T::FreudenthalTriangulation` FreudenthalTriangulation struct
- `b::Vector{Float64}` A belief in the belief space
- `bv::Vector{Float64}` The belief vertices of a triangulation `T`
- `B::Vector{Vector{Float64}}` The vector of vertices of the belief simplex of a belief `b`
- `coords::Vector{Float64}` The barycentric coordinates of the belief in the belief simplex

## Example: Interpolating a Function in Belief Space

For an example of function approximation using this package, consider the function `U(x,y) = ((x-1)^2 + (y-1)^2)^0.5`.  However, assume that we do not have access to `U(x,y)` but instead, we only know `U` for values `[x,y]` in `belief_vertices(FreudenthalTriangulation(2, m))` for some granularity `m`. Assume for this example that the maximum granularity is `m = 3`. Thus we have access to a function `U_vertices(x, y)` which returns `U(x,y)` if `[x,y]` in `belief_vertices(FreudenthalTriangulation(2, 3))` and `None` otherwise.

Using `U_vertices`, we can create a function `U_interp` which will approximate `U(x,y)` by interpolating on the Freudenthal space. To do this, we will use the `freudenthal_simplex_and_coords` function.

```julia
function U_interp(x, y)
	X = [x, y]
	T = FreudenthalTriangulation(2, 3)
	B, coords = belief_simplex(T, X)

	interp_val = 0
	for (b, coord) in zip(B, coords)
		interp_val += U_vertices(v[1], v[2]) * coord
	end
	return interp_val
end		
```
Thus `U_interp` allows us to approximate `U` at real coordinates based on the `U_vertices` function.
