using FreudenthalTriangulation
using Test

@testset "FreudenthalTriangulation.jl" begin
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

    # test freudenthal_simplex_and_coords
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

    # test to_belief and to_freudenthal
    n, m = rand(1:10, 2)
    d = sort(cat(0, rand(0:m, n-1), dims=1))
    x = [m - d[i] for i=1:n]
    b = to_belief(x, m)
    @test abs(sum(b) - 1.0) < 1e-8 && size(b, 1) == n
    @test abs(sum(to_freudenthal(b, m) - x)) < 1e-8

    # test to_freudenthal_batch and to_belief_batch
    k, m, n= rand(1:10, 3)
    B = rand(Float64, (n, k))
    B = broadcast(abs, B)
    for i = 1:k
        B[:, i] = B[:, i] ./ sum(B[:, i])
    end
    X = to_freudenthal_batch(B, m)
    for i = 1:k
        X[:, i] = to_belief(X[:, i], m)
    end
    @test abs(sum(X - B)) < 1e-8

    BB = to_belief_batch(X, m)
    @test abs(sum(BB - B) < 1e-8)
end
