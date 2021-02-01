using FreudenthalTriangulation
using Test

@testset "FreudenthalTriangulations.jl" begin
    # test freudenthal_vertices
    n, m = rand(1:10, 2)
    V = freudenthal_vertices(n, m)
    @test size(V, 1) == binomial(n + m - 1, m)
    v = rand(V)
    @test size(v, 1) == n && v[1] == m
    @test issorted(v, rev=true)

    # test freudenthal_simplex
    n = rand(5:20)
    x = rand(Int64, n)
    V = freudenthal_simplex(x)
    distances = [abs(sum(x - V[i])) for i in 1:n]
    @test size(V, 1) == n + 1
    @test issorted(distances)

    # test barycentric_coordinates
    n = rand(5:20)
    x = rand(Int64, n)
    V = freudenthal_simplex(x)
    λ = barycentric_coordinates(x, V)
    @test size(λ, 1) == n + 1
    @test sum(λ) == 1.0

    # test to_belief and to_freudenthal
    n, m = rand(1:10, 2)
    d = sort(cat(0, rand(0:m, n-1), dims=1))
    x = [m - d[i] for i=1:n]
    b = to_belief(x, m)
    @test sum(b) == 1.0 && size(b, 1) == n
    @test to_freudenthal(b, m) == x













end
