using Test, .PHysicalTDA

@testset "Crystal + LSWT + TDA smoke" begin
	cryst, sys = La2CuO4()
	@test !isnothing(cryst) && !isnothing(sys)
	
	(r2D, rpaths, Es) = PHysicalTDA.LSWT(cryst, sys)
	result2D, Qs = r2D; @test !isempty(Es) && !isempty(Qs)
	results, paths = rpaths; @test !isempty(results) && !isempty(paths)
		
	A4 = convert4D(result2D, Qs, Es); @test ndims(A4) == 4
	hkl= collapes(A4; over=:w, op=sum); @test ndims(hkl) == 3

	PD, _ = pd_sunny_intensities(result2D; maxdim=1, superlevel=false)
	@test length(PD) >= 1

	#Light determinism check on path count
	@test length(results) == length(paths)
end