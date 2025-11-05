"""
Constructs a tetragonal La₂CuO₄-like crystal and Sunny spin system.
Applies J, J′, J″, Jc exchanges, randomizes spins, and minimizes energy.

Space Group: 139 ~ Single Cu Ion @ [0,0,0] | (s=1/2, g=2)
Supercell: (6×6×3) 

Returns (cryst::Sunny.Crystal, sys::Sunny.System).
"""
function La2CuO4()
        units = Units(:meV, :angstrom)
	latvecs = lattice_vectors(3.85, 3.85, 12.25, 90, 90, 90)
	positions = [[0,0,0]]; types = ["Cu"]; sg = 139
	cryst = Crystal(latvecs, positions, sg; types)

	moments = [1=>Moment(s=1/2, g=2)]; dims = (2,2,1)
	sys = System(cryst, moments, :dipole; dims=dims)

	J = 138.3; Jp = 2; Jpp = 2; Jc = 38
	nn = Bond(1,1,[1,0,0])
	nnn =  Bond(1,1,[1,1,0])
	nnnn = Bond(1,1,[2,0,0])

	set_exchange!(sys, J-Jc/2, nn)
	set_exchange!(sys, Jp-Jc/4, nnn)
	set_exchange!(sys, Jpp, nnnn)

	randomize_spins!(sys); minimize_energy!(sys)

	return cryst, sys
end
"""
Constructs low-temperature monoclinic (CH₃)₄[NMnCl₃] crystal via Sunny.
Applies exchange coupling J and anisotropic exchange Δ (along z-axis).

This crystal is commonly referred to as TMMC. It is a 1D high-spin 
Heisenberg chain.

Space group: 14 (P2₁/b) ~ Single Mn Ion @ [1/5, 1/4, 1/5] | (s=5/2, g=2)
Supercell: (1×4×150)

Returns (cryst::Sunny.Crystal, sys::Sunny.System)
"""
function TMMC()
    units = Units(:meV, :angstrom)
    latvecs = lattice_vectors(8.811, 13.265, 6.46, 90, 99.18, 90)
    positions = [[0.2034, 0.2546, 0.2072]]; sg = 14; types = ["Mn"]
    cryst = Crystal(latvecs, positions, sg; types)
    sys = System(cryst, [1=>Moment(s=5/2, g=2)], :dipole, dims=(1,4,150))
    Δ = 0.016; J = sqrt(1/(2Δ))/5; exchange = J*Diagonal([1,1,1-Δ])
    nn = Bond(1,2,[0,0,0]); set_exchange!(sys, exchange, nn)
    randomize_spins!(sys); minimize_energy!(sys; maxiters=7500)
    return cryst, sys
end
