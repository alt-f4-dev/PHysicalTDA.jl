"""
Constructs a tetragonal La₂CuO₄-like crystal and Sunny spin system.
Applies J, J′, J″, Jc exchanges, randomizes spins, and minimizes energy.
Returns (crystal, system)

Space Group: 139 ~ Single Cu Ion @ [0,0,0] | (s=1/2, g=2)
Supercell: (6x6x3) 

Returns (cryst::Sunny.Crystal, sys::Sunny.System).
"""
function La2CuO4()
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
