using FreudenthalTriangulation
using Test

@testset "FreudenthalTriangulations.jl" begin
    n = rand(10)
    m = rand(10)
    @test size(freudenthal_vertices(n, m), 1) == factorial(n + m - 1)/(factorial(n - 1)\factorial(m))
end
