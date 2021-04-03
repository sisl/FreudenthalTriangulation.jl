module FreudenthalTriangulations

using LinearAlgebra
using SparseArrays

export FreudenthalTriangulation, vertices, simplex, belief_vertices, belief_simplex

# Construct the Freudenthal triangulation of the belief simplex

struct FreudenthalTriangulation
        n::Int
        m::Int
end

"""
    vertices(T::FreudenthalTriangulation)
Construct the list of Freudenthal vertices in an `n` dimensional space with grid resolution `m`.
The vertices are represented by a list of `n` dimensional vectors.
"""
function vertices(T::FreudenthalTriangulation)
    V = Vector{Int}[]
    v = Vector{Int}(undef, T.n)
    v[1] = T.m
    _vertices!(V, v, 2)
    return V
end

function _vertices!(V::Vector{Vector{Int64}}, v::Vector{Int64}, i::Int64)
    n = length(v)
    if i > n
        push!(V, copy(v))
        return
    end
    for k in 0 : v[i-1]
        v[i] = k
        _vertices!(V, v, i+1)
    end
end

"""
    _freudenthal_simplex(x::Vector{Float64})
Returns the list of vertices of the simplex of point `x` in the Freudenthal grid.
"""
function _freudenthal_simplex(x::Vector{Float64})
    n = length(x)
    V = Vector{Vector{Int}}(undef, n+1)
    V[1] = floor.(Int, x)
    d = x - V[1]
    p = sortperm(d, rev=true)
    for i in 2 : n+1
        V[i] = copy(V[i-1])
        V[i][p[i-1]] += 1
    end
    return V
end

"""
    _barycentric_coordinates(x::Vector{Float64}, V::Vector{Vector{Int64}})
Given a point `x` and its simplex `V` in the Freudenthal grid, returns the barycentric coordinates
of `x` in the grid. `V` must be in the same order as provided by the output of `freudenthal_simplex`
"""
function _barycentric_coordinates(x::Vector{Float64}, V::Vector{Vector{Int64}})
    d = x - V[1]
    p = sortperm(d, rev=true)
    n = length(x)
    λ = Vector{Float64}(undef, n+1)
    λ[n+1] = d[p[n]]
    for i in n:-1:2
        λ[i] = d[p[i-1]] - d[p[i]]
    end
    λ[1] = 1.0 - sum(λ[2:end])
    return λ
end

"""
    simplex(T::FreudenthalTriangulation, x::Vector{Float64})
Given a point `x`, returns the simplex of the point `x` and the barycentric coordinates of `x` in the grid.
"""
function simplex(T::FreudenthalTriangulation, x::Vector{Float64})
    V = _freudenthal_simplex(x)
    return V, _barycentric_coordinates(x, V)
end

"""
    _to_belief(x)
Transform a point `x` in the Freudenthal space to a point in the belief space.
"""
_to_belief(x) = (push!(x[1:end-1] - x[2:end], x[end]))./x[1]

"""
    _to_freudenthal(b, m::Int64)
Transform a point `b` in the belief space to a point in the Freudenthal space.
`m` is the resolution of the Freudenthal grid.
"""
_to_freudenthal(b, m::Int64) = [sum(b[k] for k in i : length(b))*m for i in 1 : length(b)]

"""
    belief_vertices(T::FreudenthalTriangulation)
    Converts the list of Freudenthal vertices with dimension `T.n` and granularity `T.m` into
    belief space and returns the list
"""
belief_vertices(T::FreudenthalTriangulation) = _to_belief.(vertices(T))

"""
    belief_simplex(T::FreudenthalTriangulation, b)
    Constructs the belief simplex surrounding a belief `b` and the barycentric coordinates
    corresponding to the location of the belief within the simplex
"""
function belief_simplex(T::FreudenthalTriangulation, b)
    x = _to_freudenthal(b, T.m)
    V, λ = simplex(T, x)
    B = _to_belief.(V)
    valid =  λ .> sqrt(eps())
    return B[valid], λ[valid]
end

end # module
