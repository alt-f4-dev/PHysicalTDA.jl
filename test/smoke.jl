using Test, PHysicalTDA

@testset "Crystal + LSWT + TDA smoke" begin
    cryst, sys = La2CuO4()
    @test !isnothing(cryst) && !isnothing(sys)
    
    Nq = 500
    h = [[-1,0,0], [0,0,0], [1,0,0]]
    k = [[0,-1,0], [0,0,0], [0,1,0]]
    l = [[0,0,-1], [0,0,0], [0,0,1]]
    H = Sunny.q_space_path(cryst, h, Nq)
    K = Sunny.q_space_path(cryst, k, Nq)
    L = Sunny.q_space_path(cryst, l, Nq)
    paths = [H,K,L]
    LCO = PHysicalTDA.LSWT(cryst, sys, paths)
    qpaths = LCO.qpaths; ωaxis = LCO.energies
    qpathIntensities = LCO.intensities
    @test !isempty(qpaths) && !isempty(ωaxis)
    @test !isempty(qpathIntensities)
     
    PD = pd_sunny_intensities(qpathIntensities[1]; 
                              maxdim=1, superlevel=false)
    @test length(PD) ≥ 1

    #Light determinism check on path count
    @test length(qpathIntensities) == length(qpaths)
    qpathIntensity1 = qpathIntensities[1].data; qpath1 = qpaths[1].qs
    @test size(qpathIntensity1) == (length(qpath1), length(ωaxis))
    
    #To Add: 
    #       - TopoViz.jl tests
    #       - Fingerprint tests
end
