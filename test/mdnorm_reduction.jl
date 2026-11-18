using LinearAlgebra

function reduction_histogram(centres, D, N; D2=copy(N), names=nothing, cell=nothing)
    coords = [collect(Float64, axis) for axis in centres]
    edges = map(coords) do x
        length(x) == 1 && return [x[1]-0.5, x[1]+0.5]
        [(3x[1]-x[2])/2; (x[1:end-1] .+ x[2:end]) ./ 2; (3x[end]-x[end-1])/2]
    end
    names = names === nothing ? collect(("H", "K", "DeltaE", "L")[1:length(centres)]) : collect(names)
    units = [name == "DeltaE" ? "DeltaE" : "r.l.u." for name in names]
    frames = [name == "DeltaE" ? "General Frame" : "HKL" for name in names]
    meta = cell === nothing ? (; runs=[]) : (; runs=[(; name="experiment0", lattice=cell)])
    return NormalizedHistogram(D, D2, N, edges, coords, names, units, frames, meta)
end

@testset "Exposure-weighted MDNorm reductions" begin
    D = [2.0 10; 8 20; 900 21]
    N = [1.0 2; 2 4; 0 3]
    D2 = [1.0 2; 4 16; 500 9]
    h = reduction_histogram(([0.0,1,2],[0.5,1.5]), D, N; D2, names=("H","DeltaE"))
    b = rebin(h; factors=(2,1))
    @test b.D == [10.0 30.0]
    @test b.N == [3.0 6.0]
    @test b.D2 == [5.0 18.0]
    @test b.edges[1] == [-0.5,1.5]
    @test intensity(b) == [10/3 5]
    @test uncertainty(b) ≈ [sqrt(5)/3 sqrt(18)/6]
    keep = rebin(h; factors=(2,1), remainder=:keep)
    @test keep.N == [3.0 6.0; 0 3]
    @test keep.D[2,1] == keep.D2[2,1] == 0
    @test isnan(intensity(keep)[2,1])
    @test isnan(uncertainty(keep)[2,1])
    @test rebin(h; factors=(2,1), remainder=:merge).D == [10.0 51.0]
    @test rebin_to_resolution(h; factors=(2,1)).D == b.D
    @test rebin_to_resolution(h; resolution_q=2, resolution_ω=1).D == b.D
    x = integrate(h; dims=2, range=(0.5,1.5), weight=i -> i[2]==1 ? 2.0 : 0.5)
    @test vec(x.D) == [9.0,26.0,10.5]
    @test vec(x.N) == [3.0,6.0,1.5]
    @test vec(x.D2) == [4.5,20.0,2.25]
    @test x.centres[2] == [1.0]
    @test integrate(h; dims=(1,2), range=(1:2,2:2)).D == fill(30.0,1,1)
    @test_throws ArgumentError integrate(h; dims=1, range=(4.0,5.0))
    @test_throws ArgumentError integrate(h; dims=2, weight=-1)
    @test_throws DimensionMismatch rebin(h; factors=(3,))
    h_noerr = reduction_histogram(h.centres, D, N; D2=nothing, names=h.names)
    @test uncertainty(rebin(h_noerr; factors=2)) === nothing
    @test weighted_average(h,[1.0 2; 4 4; NaN 9]; dims=1) ≈ [3.0 47/9]
    @test isnan(only(weighted_average(h,NaN; dims=(1,2))))
    target = reduction_histogram(([0.5,1.1,2.6],[-0.2,1.0,1.8]), zeros(3,3), ones(3,3); names=h.names)
    f = reshape([3.0,5.0,NaN],3,1)
    mapped = map_nearest_bins(f,h,target)
    @test mapped.values[1,:] == fill(3.0,3)
    @test mapped.values[2,:] == fill(5.0,3)
    @test all(isnan,mapped.values[3,:])
    @test !any(mapped.valid[3,:])
    @test_throws DimensionMismatch map_nearest_bins(f,h,target; dims=(1,2))
end

@testset "Bragg indexing correction" begin
    q = collect(-1.0:0.05:1.0)
    M = [1.0 0.1; -0.1 1.0]
    t = [0.05,-0.05]
    P = [(-0.5,-0.5),(0.5,-0.5),(-0.5,0.5),(0.5,0.5)]
    N = [1.0 + mod(i+2j,4) for i in eachindex(q),j in eachindex(q),k in 1:2]
    D = copy(N)
    for p in P
        x,y = M*collect(p)+t
        i,j = argmin(abs.(q .- x)), argmin(abs.(q .- y))
        D[i,j,1] = 20N[i,j,1]
    end
    D[:,:,2] .= 200 .* N[:,:,2]
    h = reduction_histogram((q,q,[0.0,1.0]),D,N; D2=0.01 .* N.^2)
    peaks = find_bragg_peaks(h; positions=P, energy_window=(-0.1,0.1), window=0.16)
    @test length(peaks) == 4
    @test all(r -> r.contrast>=5,peaks)
    fit = fit_index_correction(peaks)
    @test fit.M ≈ M atol=1e-14
    @test fit.t ≈ t atol=1e-14
    @test fit.rms < 1e-14 < fit.rms_uncorrected
    @test collect(corrected_indices(h,fit.M,fit.t)[1,end])[1:2] ≈ M \ ([first(q),last(q)]-t)
    @test corrected_indices(h,Matrix{Float64}(I,2,2),[0.0,0.0])[4,4] == (q[4],q[4],0.0)
    line = [(; P=(Float64(x),0.0), observed=(Float64(x),0.0), contrast=1.0) for x in 1:3]
    @test_throws ArgumentError fit_index_correction(line)
    @test_throws ArgumentError corrected_indices(h,zeros(2,2),t)
    @test isempty(find_bragg_peaks(h; positions=[(0.0,0.0)],energy_window=1:1))
    M3 = [1.0 0.04 0.02; 0.0 1.03 -0.02; 0.01 0.0 0.98]
    t3 = [0.04,-0.07,0.01]
    p3 = [(0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0)]
    r3 = [(;P=p,observed=Tuple(M3*collect(p)+t3),contrast=1.0) for p in p3]
    f3 = fit_index_correction(r3; weights=:uniform)
    @test f3.M ≈ M3
    @test f3.t ≈ t3
end

@testset "Reciprocal cuts and shared source covariance" begin
    q = [-0.375,-0.125,0.125,0.375]
    shape = (4,4,3,1)
    D,D2,N = zeros(shape),zeros(shape),zeros(shape)
    for i in 1:4,e in 1:3
        N[i,i,e,1] = i
        D[i,i,e,1] = i * ((i in (2,3) ? 2 : 3)==e ? 9 : 1)
        D2[i,i,e,1] = i^2
    end
    h = reduction_histogram((q,q,[1.0,2.0,3.0],[0.25]),D,N; D2)
    eye = Matrix{Float64}(I,2,2)
    group = (eye,-eye)
    path = [(0.0,0.0),(0.5,0.5)]
    args = (;width=0.01,step=sqrt(2)/4)
    unfolded = cut_path(h,path;args...)
    folded = cut_path(h,path;args...,fold=group)
    @test unfolded.values == folded.values == [1.0 9 1; 1 1 9]
    @test unfolded.weight == [3.0 3 3; 4 4 4]
    @test folded.weight == fill(5.0,2,3)
    @test isempty(folded.covariance)
    @test folded.σ ≈ repeat([sqrt(13)/5,sqrt(17)/5],1,3)
    @test !any(cut_path(h,path;args...,indices=corrected_indices(h,eye,[0.1,0.0])).valid)
    mask = falses(shape)
    mask[3,3,:,:] .= true
    @test !any(cut_path(h,path;args...,exclude=mask).valid[1,:])
    @test !any(cut_path(h,path;args...,exclude=q -> q[1]>0.25).valid[2,:])
    @test symmetry_images((0.5,0.0),group) == [(-0.5,0.0)]
    @test Set(symmetry_images((1.125,0.125),group)) == Set([(0.125,0.125),(-0.125,-0.125)])
    swap = [0.0 1; 1 0]
    @test Set(symmetry_images((0.1,0.2),(eye,swap);translations=false)) == Set([(0.1,0.2),(0.2,0.1)])
    eye3 = Matrix{Float64}(I,3,3)
    @test symmetry_images((0.0,0.25,1.5),(eye3,)) == [(0.0,0.25,-0.5)]
    single = reduction_histogram(([0.25],[0.25],[1.0]),fill(6.0,1,1,1),fill(3.0,1,1,1);D2=fill(4.0,1,1,1))
    grid = [-0.5,0.0,0.5]
    zone = fold_to_zone(single,group;grid=(grid,grid))
    @test zone.N[:,:,1] == [3.0 0; 0 3]
    @test zone.D[:,:,1] == [6.0 0; 0 6]
    @test zone.meta.numerator_covariance == [(1,4,4.0)]
    combined = integrate(zone;dims=(1,2))
    @test only(combined.D2) == 16
    @test only(uncertainty(combined)) ≈ 2/3
    @test rebin(zone; factors=(2,2,1)).D2 == combined.D2
    weighted = integrate(zone; dims=(1,2), weight=i -> i[1]==1 ? 2.0 : 3.0)
    @test only(weighted.D2) == 100
    @test only(uncertainty(weighted)) ≈ 2/3
    zone3 = fold_to_zone(single,(eye3,);grid=(grid,grid,[-0.5,0.5]))
    @test size(zone3.N) == (2,2,1,1)
    @test zone3.N[2,2,1,1] == 3
    across = cut_path(single,[(-0.5,-0.5),(0.5,0.5)];args...,fold=group)
    @test length(across.covariance) == 1
    @test only(across.covariance)[3] ≈ 4/9
    @test only(fold_to_zone(zone,eye;grid=([-0.5,0.5],[-0.5,0.5])).D2) == 16
    collapsed = fold_to_zone(single,(eye,eye,-eye);grid=([-0.5,0.5],[-0.5,0.5]))
    @test only(collapsed.N) == 3
    @test only(collapsed.D2) == 4
    @test isempty(collapsed.meta.numerator_covariance)
    @test_throws ArgumentError fold_to_zone(single,group;grid=([-0.5,0.0,0.0,0.5],grid))
end

@testset "Magnetic conversion and L averaging" begin
    z = (0,0,1)
    @test polarization_factor(z,z) == 1
    @test polarization_factor((2,0,0),z) == 1/2
    @test polarization_factor(z,z;channel=:longitudinal) == 0
    @test polarization_factor((1,0,0),z;channel=:longitudinal) == 1
    @test polarization_factor((0,0,0),z;channel=:isotropic) == 2/3
    @test isnan(polarization_factor((0,0,0),z))
    N = reshape([1.0,3.0],1,1,1,2)
    coords = ([1.0],[0.0],[2.0],[0.0,1.0])
    cell = (;a=2π,b=2π,c=2π,alpha=90.0,beta=90.0,gamma=90.0)
    expected = reshape([1/2,3/4],1,1,1,2)
    D = @. N*(2+6expected)
    h = reduction_histogram(coords,D,N;cell,D2=fill(4.0,size(N)))
    g = magnetic_conversion(h;ion=q -> 1.0)
    @test g ≈ expected
    avg = magnetic_conversion(h;ion=q -> 1.0,average_over=4)
    @test size(avg) == (1,1,1,1)
    @test only(avg) ≈ (1/2+3*3/4)/4
    @test avg ≈ weighted_average(h,g;dims=4)
    @test magnetic_conversion(h;ion=q -> 1.0,channel=:isotropic) ≈ fill(2/3,size(N))
    @test magnetic_conversion(h;ion=:Mn2)[1] ≈ magnetic_form_factor(1.0,MN2_J0)^2/2
    @test magnetic_conversion(h;ion=q -> 1.0,channel=:longitudinal) ≈ 2 .* (1 .- g)
    background = (q,ω) -> (2.0,0.5)
    S,σ = to_smag(h,g;background)
    @test all(x -> x≈6,S)
    @test σ ≈ sqrt.(h.D2 .+ (0.5 .* N).^2) ./ (N .* g)
    path = [(0.5,0.0),(1.5,0.0)]
    cut = cut_path(h,path;width=0,step=1,conv=g,background)
    @test only(cut.values) ≈ 6
    @test only(cut.weight) ≈ 2.75
    @test only(cut.σ) ≈ sqrt(8+0.25+2.25)/2.75
    @test isnan(to_smag(h,0.0)[1][1])
    @test_throws ArgumentError to_smag(h,-1.0)
    @test_throws DimensionMismatch to_smag(h,[1.0])
    @test cut_path(h,path;width=0,step=1,conv=NaN).valid == falses(1,1)
    skew = merge(cell,(;gamma=60.0))
    hs = reduction_histogram(coords,D,N;cell=skew)
    @test magnetic_conversion(hs;ion=identity,channel=:isotropic)[1] ≈ (2/3)/sind(60)^2
    @test to_smag(reduction_histogram(coords,D,N;D2=nothing),g)[2] === nothing
end

@testset "MDNorm disk-to-cut A/B workflow" begin
    mktempdir() do dir
        grid_B = ([0.0,1,2],[0.0,1,2],[0.0,1,2],[0.0,1,2])
        N_B = [1.0+l for i in 1:2,j in 1:2,e in 1:2,l in 1:2]
        pb = write_synthetic_mdnorm(mkpath(joinpath(dir,"B")),copy(N_B),copy(N_B),N_B;grid=grid_B)
        B = load_mdnorm(pb.path_D,pb.path_N)
        basis = Matrix{Float64}(I,3,3)
        g_B = magnetic_conversion(B;ion=q -> 1.0,average_over=4,reciprocal_basis=basis)
        grid_A = ([0.1,0.7],grid_B[2],grid_B[3])
        N_A = reshape([1.0,5.0,2.0,4.0],1,2,2)
        ref = map([0.5,1.5]) do k
            f(l) = (1+l^2/(0.5^2+k^2+l^2))/2
            (2f(0.5)+3f(1.5))/5
        end
        S = reshape([9.0,4.0,3.0,8.0],1,2,2)
        exposure = N_A .* reshape(ref,1,2,1)
        D_A = @. S*exposure+N_A*1.25
        pa = write_synthetic_mdnorm(dir,D_A,copy(N_A),N_A;grid=grid_A)
        A = load_mdnorm(pa.path_D,pa.path_N;path_I=pa.path_I)
        mapped = map_nearest_bins(g_B,B,A)
        @test all(mapped.valid)
        @test mapped.values ≈ repeat(reshape(ref,1,2,1),1,1,2)
        hkl = corrected_indices(A,basis[1:2,1:2],[-0.1,0.0])
        cut = cut_path(A,[(0.5,0.0),(0.5,2.0)];width=0,step=1,indices=hkl,conv=mapped.values,background=1.25)
        @test cut.values ≈ reshape(S,2,2)
        @test cut.weight ≈ reshape(exposure,2,2)
        @test cut.σ ≈ reshape(sqrt.(N_A) ./ exposure,2,2)
    end
end
