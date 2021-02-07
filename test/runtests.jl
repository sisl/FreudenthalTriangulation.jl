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
    x = rand(Float64, n)
    V = freudenthal_simplex(x)
    @test size(V, 1) == n + 1
    x = rand(Int64, n)
    x = convert(Vector{Float64}, x)
    V = freudenthal_simplex(x)
    distances = [abs(sum(x - V[i])) for i in 1:n]
    @test issorted(distances)

    # test barycentric_coordinates
    n = rand(5:20)
    x = rand(Float64, n)
    V = freudenthal_simplex(x)
    λ = barycentric_coordinates(x, V)
    @test size(λ, 1) == n + 1
    @test abs(sum(λ) - 1.0) < 1e-8

    # test freudenthal_simplex_and_coords!
    n = rand(5:20)
    x = rand(Float64, n)
    V = freudenthal_simplex(x)
    λ = barycentric_coordinates(x, V)
    V2, λ2 = freudenthal_simplex_and_coords(x)
    i = rand(1:n)
    v = V[i]
    v2 = V2[i]
    @test abs(sum(v - v2)) < 1e-8
    @test abs(sum(λ - λ2)) < 1e-8

    # test freudenthal_simplex_and_coords!
    n = rand(5:20)
    x = rand(Float64, n)
    V = freudenthal_simplex(x)
    λ = barycentric_coordinates(x, V)
    V2 = [rand(Int64, n) for i=1:n+1]
    λ2 = rand(Float64, n+1)
    V2, λ2 = freudenthal_simplex_and_coords!(x, V2, λ2)
    i = rand(1:n)
    v = V[i]
    v2 = V2[i]
    @test abs(sum(v - v2)) < 1e-8
    @test abs(sum(λ - λ2)) < 1e-8
end
