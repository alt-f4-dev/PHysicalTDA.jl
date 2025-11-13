using Test, PHysicalTDA

@testset "Crystal + LSWT + TDA smoke" begin
	cryst, sys = La2CuO4()
	@test !isnothing(cryst) && !isnothing(sys)
	
	LCO = PHysicalTDA.LSWT(cryst, sys)
        qpaths = LCO.qpaths; ωaxis = LCO.energies
        qpathIntensities = LCO.intensities
	@test !isempty(qpaths) && !isempty(ωaxis)
	@test !isempty(qpathIntensities)
        
	PD = pd_sunny_intensities(result2D; maxdim=1, superlevel=false)
	@test length(PD) ≥ 1

	#Light determinism check on path count
	@test length(qpathIntensities) == length(qpaths)
        
        #To Add: 
        #       - TopoViz.jl tests
        #       - Fingerprint tests
end
