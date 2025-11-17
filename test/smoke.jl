using Test, Sunny, PHysicalTDA

@testset "Crystal + LSWT + TDA smoke" begin
    #------------------#
    # Crystal Database #
    #------------------#
    cryst, sys = La2CuO4()
    @test !isnothing(cryst) && !isnothing(sys)
    
    #--------------#
    # LSWTtools.jl #
    #--------------#
    Nq = 500
    h = [[-1,0,0], [0,0,0], [1,0,0]]
    k = [[0,-1,0], [0,0,0], [0,1,0]]
    l = [[0,0,-1], [0,0,0], [0,0,1]]
    H = q_space_path(cryst, h, Nq)
    K = q_space_path(cryst, k, Nq)
    L = q_space_path(cryst, l, Nq)
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
    
    #----------------------------#
    # Filtration & Normalization #
    #----------------------------#
    A = [0.0 1.0;
         0.0 2.0]
    
    #Sublevel: (normalize, threshold) ~ (false, nothing) ⟹   Z ∈ [min(A), max(A)] & θ = 1.0
    PD_sub_full_unnorm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=false, threshold=nothing)
    PD_sub_filtered_unnorm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=false, threshold=1.0)
    @test length(PD_sub_full_unnorm[1]) == length(PD_sub_filtered_unnorm[1])
    b1, d1 = PD_sub_full_unnorm[1].intervals[1]
    b2, d2 = PD_sub_filtered_unnorm[1].intervals[1]
    @test (b1 == b2) && (d1 == d2)
    
    #Superlevel: (normalize, threshold) ~ (false, nothing) ⟹   Z ∈ [min(A), max(A)] & θ = 0.0
    PD_sup_full_unnorm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=false, threshold=nothing)
    PD_sup_filtered_unnorm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=false, threshold=0.0)
    @test length(PD_sup_full_unnorm[1]) == length(PD_sub_filtered_unnorm[1])
    b1, d1 = PD_sup_full_unnorm[1].intervals[1]
    b2, d2 = PD_sup_filtered_unnorm[1].intervals[1]
    @test (b1 == b2) && (d1 == d2)

    #Sublevel: (normalize, threshold) ~ (true, nothing) ⟹   Z ∈ [0, 1] & θ = 1.0
    PD_sub_full_norm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=true, threshold=nothing)
    PD_sub_filtered_norm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=true, threshold=1.0)
    @test length(PD_sub_full_norm[1]) == length(PD_sub_filtered_norm[1])
    b1, d1 = PD_sub_full_norm[1].intervals[1]
    b2, d2 = PD_sub_filtered_norm[1].intervals[1]
    @test (b1 == b2) && (d1 == d2)

    #Superlevel: (normalize, threshold) ~ (true, nothing) ⟹   Z ∈ [0, 1] & θ = 0.0
    PD_sup_full_norm = pd_array_intensities(A; maxdim=0, superlevel=true, normalize=true, threshold=nothing)
    PD_sup_filtered_norm = pd_array_intensities(A; maxdim=0, superlevel=false, normalize=true, threshold=0.0)
    @test length(PD_sup_full_norm[1]) == length(PD_sup_filtered_norm[1])
    b1, d1 = PD_sup_full_norm[1].intervals[1]
    b2, d2 = PD_sup_filtered_norm[1].intervals[1]
    @test (b1 == b2) && (d1 == d2) 

    #----------------------------#
    # To Add:                    #
    #        - TopoViz.jl tests  #
    #        - Fingerprint tests #
    #----------------------------#
end
