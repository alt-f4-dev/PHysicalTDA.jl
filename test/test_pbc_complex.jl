#--------------------------------------#
#  PeriodicCubical Verification Tests  #
#--------------------------------------#----------------------------------------#
#                                                                               #
# Ground-truth homology:                                                        #
#   S¹  (circle)      H₀=ℤ,  H₁=ℤ                                               #
#   T²  (2-torus)     H₀=ℤ,  H₁=ℤ²,  H₂=ℤ                                       #
#   S¹×[0,1]          H₀=ℤ,  H₁=ℤ         (cylinder)                            #
#   [0,1]²            H₀=ℤ                 (contractible)                       #
#                                                                               #
# Constant filtration: all values equal → essential classes = Betti numbers.   #
#-------------------------------------------------------------------------------#
@testset "PeriodicCubical verification" begin

    function count_essential(dgms)
        map(dgms) do d
            count(p -> isinf(PersistenceDiagrams.death(p)), d)
        end
    end

    function count_finite(dgms)
        map(dgms) do d
            count(p -> isfinite(PersistenceDiagrams.death(p)), d)
        end
    end

    # Test 1: Aperiodic square (contractible)
    @testset "Aperiodic square (contractible)" begin
        img  = ones(Float64, 10, 10)
        cf   = PeriodicCubical(img; periodic_dims=(false, false))
        dgms = ripserer(cf; dim_max=2)
        ess  = count_essential(dgms)
        @test ess[1] == 1
        @test ess[2] == 0
        @test ess[3] == 0
    end

    # Test 2: Circle S¹
    @testset "Circle S¹" begin
        img  = ones(Float64, 1, 20)
        cf   = PeriodicCubical(img; periodic_dims=(false, true))
        dgms = ripserer(cf; dim_max=1)
        ess  = count_essential(dgms)
        @test ess[1] == 1
        @test ess[2] == 1
    end

    # Test 3: Cylinder S¹×[0,1] — physically relevant BZ × energy strip
    @testset "Cylinder S¹×[0,1] (periodic q, aperiodic ω)" begin
        img  = ones(Float64, 20, 10)
        cf   = PeriodicCubical(img; periodic_dims=(true, false))
        dgms = ripserer(cf; dim_max=2)
        ess  = count_essential(dgms)
        @test ess[1] == 1
        @test ess[2] == 1
        @test ess[3] == 0
    end

    # Test 4: 2-Torus T²
    @testset "2-Torus T²" begin
        img  = ones(Float64, 20, 20)
        cf   = PeriodicCubical(img; periodic_dims=(true, true))
        dgms = ripserer(cf; dim_max=2)
        ess  = count_essential(dgms)
        @test ess[1] == 1
        @test ess[2] == 2
        @test ess[3] == 1
    end

    # Test 5: Tiling and cubemap shape correctness
    @testset "Cubemap shape" begin
        img = rand(Float64, 8, 12)

        t_none = PHysicalTDA._tile_periodic(img, (false, false))
        @test size(t_none) == (8, 12)
        @test size(PHysicalTDA._periodic_cubemap(img, (false, false))) == (15, 23)

        t1 = PHysicalTDA._tile_periodic(img, (true, false))
        @test size(t1) == (9, 12)
        @test size(PHysicalTDA._periodic_cubemap(img, (true, false))) == (16, 23)

        t2 = PHysicalTDA._tile_periodic(img, (false, true))
        @test size(t2) == (8, 13)
        @test size(PHysicalTDA._periodic_cubemap(img, (false, true))) == (15, 24)

        t12 = PHysicalTDA._tile_periodic(img, (true, true))
        @test size(t12) == (9, 13)
        @test size(PHysicalTDA._periodic_cubemap(img, (true, true))) == (16, 24)

        @test t1[9, :]  == t1[1, :]
        @test t2[:, 13] == t2[:, 1]
    end

    # Test 6: Consistency with Cubical when all aperiodic
    @testset "PeriodicCubical == Cubical when all aperiodic" begin
        img      = rand(Float64, 10, 15)
        dgms_ref = ripserer(Cubical(img); dim_max=1)
        cf_per   = PeriodicCubical(img; periodic_dims=(false, false))
        dgms_per = ripserer(cf_per; dim_max=1)
        for d in 1:2
            @test length(dgms_ref[d]) == length(dgms_per[d])
        end
    end

    # Test 7: pd_array_intensities integration
    @testset "pd_array_intensities with periodic_dims" begin
        img   = rand(Float64, 20, 50)
        dgms1 = pd_array_intensities(img; maxdim=1, periodic_dims=nothing)
        dgms2 = pd_array_intensities(img; maxdim=1, periodic_dims=(false, false))
        dgms3 = pd_array_intensities(img; maxdim=1, periodic_dims=(true, false))
        @test length(dgms1) == 2
        @test length(dgms2) == 2
        @test length(dgms3) == 2
        ess_aperiodic  = count(p -> isinf(PersistenceDiagrams.death(p)), dgms2[2])
        ess_periodic_q = count(p -> isinf(PersistenceDiagrams.death(p)), dgms3[2])
        @test ess_periodic_q ≥ ess_aperiodic + 1
    end

    # Test 8: Superlevel filtration with PBCs
    @testset "Superlevel + periodic" begin
        img  = rand(Float64, 15, 15)
        dgms = pd_array_intensities(
            img; maxdim=1, superlevel=true, periodic_dims=(true, true)
        )
        @test length(dgms) == 2
        ess = count(p -> isinf(PersistenceDiagrams.death(p)), dgms[2])
        @test ess ≥ 2
    end
end
