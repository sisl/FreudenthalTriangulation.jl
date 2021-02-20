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

## Example

For an example of function approximation using this package, consider the function `f(x,y) = ((x-1)^2 + (y-1)^2)^0.5`.  However, assume that we do not have access to `f(x,y)` but instead, we only know `f` at integer values of `x` and `y` so we have access to a function `f_integer(x, y)` which returns `f(x,y)` if `x` and `y` are integers and `None` otherwise.

Using `f_integer`, we can create a function `f_interp` which will approximate `f(x,y)` by interpolating on the Freudenthal space. To do this, we will use the `freudenthal_simplex_and_coords` function.

```julia
function f_interp(x, y)
	X = [x, y]
	V, coords = freudenthal_simplex_and_coords(X)

	interp_val = 0
	for (v, coord) in zip(V, coords)
		interp_val += f_integer(v[1], v[2]) * coord
	end
	return interp_val
end		
```
Thus `f_interp` allows us to approximate `f` at real coordinates based on the `f_integer` function.
