"""
Computes LSWT S(q,ω) on along several q-space paths (uses `Sunny.jl`). 
Build broadened 2D intensity and multiple 1D path scans in parallel.
The default kernel is Gaussian, if not provided.

Input: `cryst::Sunny.Crystal, sys::Sunny.System, paths::Vector{Sunny.QPath}`
Output: `(intensities, qpaths, energies)::NamedTuple`
"""
function LSWT(cryst, sys, paths; energies = range(0,350,500), kern = gaussian(fwhm=35))
    results = Vector{Sunny.Intensities{Float64, Sunny.QPath, 2}}(undef, length(paths))
    @threads for i in eachindex(paths)
	    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
	    results[i] = intensities(swt, paths[i]; energies, kernel=kern)
    end
    return (intensities = results, qpaths = paths, energies = energies)
end
"""
Reshape Sunny intensities from (ω, q) or (q, ω) into a 4D cube (h,k,ℓ,ω).
Validates sizes given Q-grid shape, then reshapes/permutedims accordingly.
Error thrown if flattened data size does not match nω × nQ.
Returns a 4D Array with axis order (h,k,ℓ,ω).
"""
function convert4D(result2D, Qs, Es)#Output ~ (h,k,l,w)
	A = Array(result2D.data) #Dimensions ~ (w,q)
	nω = length(Es); nq = size(Qs); nQ = prod(nq)

	@assert length(A) == nω * nQ "Data size mismatch: got $(size(A)), expected $(nω) × $(nQ)"
	if size(A,1) == nω && size(A,2) == nQ
		return permutedims(reshape(A, (nω, nq...)), (2,3,4,1))
	elseif size(A,1) == nQ && size(A,2) == nω
		return reshape(A, (nq...,nω))
	else 
		error("Unexpected Dimension in Convert4D(): $(size(A))")
	end
end
"""
Collapse (reduce) a 4D array over one or more axes with specified operation (default operation is :sum).
Axes may be symbols (:h,:k,:ℓ,:ω) [Julia Syntax] or integer indices [Python Syntax].
Reduces highest axes first to keep indices stable. Returns collapsed (reduced) array.
 
Note, the indexing convention in Julia is such that (JuliaArr[1] ~ PythonArr[0]).
"""
const AX = (; h=1, k=2, ℓ=3, l=3, ω=4, w=4)
function collapse(A; over=:ω, op=sum)
	axes = over isa Tuple ? over : (over,)
    	idxs = sort!(map(x -> x isa Symbol ? AX[x] : Int(x), collect(axes)))
    	for ax in Iterators.reverse(idxs)       # reduce highest axis first
        	A = dropdims(op(A; dims=ax); dims=ax)
    	end
	return A
end
