using Test, PHysicalTDA
include("smoke.jl")

@testset "Periodic Cubical Complex" begin
    include("test_pbc_complex.jl")
end
