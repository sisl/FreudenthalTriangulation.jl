using FreudenthalTriangulations
using Test

@testset "FreudenthalTriangulations.jl" begin
    n, m = rand(1:10, 2)
    T = FreudenthalTriangulation(n, m)

    # test vertices
    V = vertices(T)
    @test size(V, 1) == binomial(n + m - 1, m)
    v = rand(V)
    @test size(v, 1) == n && v[1] == m
    @test issorted(v, rev=true)

    # test simplex
    x = rand(Float64, T.n)
    V, λ = simplex(T, x)
    @test size(V, 1) == n + 1
    @test size(λ, 1) == n + 1
    @test abs(sum(λ) - 1.0) < 1e-8
    x = rand(Int64, T.n)
    x = convert(Vector{Float64}, x)
    V, λ = simplex(T, x)
    distances = [abs(sum(x - V[i])) for i in 1:n]
    @test issorted(distances)

    # test belief_vertices
    B = belief_vertices(T)
    b = rand(B)
    @test abs(sum(b) - 1.0) < 1e-8

    # test belief_simplex
    b = rand(Float64, T.n)
    b = b ./ sum(b)
    B, λ = belief_simplex(T, b)
    bv = rand(B)
    @test abs(sum(bv) - 1.0) < 1e-8
    @test abs(sum(λ) - 1.0) < 1e-8
end
