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

Finally, to convert a batch of points in freudenthal space `X` into belief space with granularity `m`, use
```julia
X = to_belief_batch(X, m)
```
For these functions the requirements are
- `m::Int64` Granularity of the Freudenthal triangulation
- `B::AbstractArray` The columns are the belief state points
- `X::AbstractArray` The columns are the points in Freudenthal space

## Example: Interpolating a Function in Freudenthal Space

For an example of function approximation using this package, consider the function `f(x,y) = ((x-1)^2 + (y-1)^2)^0.5`.  However, assume that we do not have access to `f(x,y)` but instead, we only know `f` at integer values of `x` and `y` so we have access to a function `f_vertices(x, y)` which returns `f(x,y)` if `x` and `y` are integers and `None` otherwise.

Using `f_vertices`, we can create a function `f_interp` which will approximate `f(x,y)` by interpolating on the Freudenthal space. To do this, we will use the `freudenthal_simplex_and_coords` function.

```julia
function f_interp(x, y)
	X = [x, y]
	V, coords = freudenthal_simplex_and_coords(X)

	interp_val = 0
	for (v, coord) in zip(V, coords)
		interp_val += f_vertices(v[1], v[2]) * coord
	end
	return interp_val
end		
```
Thus `f_interp` allows us to approximate `f` at real coordinates based on the `f_vertices` function.

## Example: Interpolating a Value Function over a Belief Space

Consider a value function over a belief space `U(b)`. We want to interpolate this function over a Freudenthal space of dimension `n = 3` and granularity `m = 3`. We can find the belief space vertices that we need to complete the interpolation.
```julia
belief_space_vert = to_belief.(freudenthal_simplex(n, m), m)
```
Then assume that for each vertex in `belief_space_vert`, we know the value function at that point in belief space. Assume that we have a function `U_vertices` which takes in a vertex in a vector of size `n`, `bv`, and returns `None` if the vector is not in `belief_space_vert` and `U(bv)` if `v` is in `belief_space_vert`. Thus we only need access to the `U_vertices` function instead of needing the value function over the entire belief space.

Now we want to interpolate the value function to any point `b`.
```julia
function U_interp(b)
 x = to_freudenthal(b, m)
 V, coords = freudenthal_simplex_and_coords(x)

 interp_val = 0
 for (v, coord) in zip(V, coords):
	 interp_val += U_vertices(to_belief(v, m)) * coord
 end
 return interp_val
end		
```
Thus `U_interp` allows us to approximate the value function over the entire belief space. 
