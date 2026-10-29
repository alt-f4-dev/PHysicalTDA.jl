using Test
using PHysicalTDA

@testset "Plotting extension is optional" begin
    # Makie is a weak dependency: loading PHysicalTDA alone must not load it.
    if !any(m -> nameof(m) === :Makie, values(Base.loaded_modules))
        @test PHysicalTDA._viz_ext() === nothing
        @test_logs (:info,) (:warn,) match_mode=:any @test_throws ErrorException enable_visuals(backend=:cairo)
    else
        @test PHysicalTDA._viz_ext() !== nothing
    end
    @test_throws ErrorException enable_visuals(backend=:bogus)
end
