"""
Computes LSWT S(q,ω) on a coarse 3D Q-grid and along several q-space paths.
Build Gaussian-broadened 2D intensity and multiple 1D path scans in parallel.
Returns ([results2D, Qs], [result_along_paths, paths], energies)

Assumes input is a Sunny.System corresponding to the provided Sunny.Crystal type.
"""
function LSWT(cryst, sys; energies = range(0,350,500), kern = gaussian(fwhm=35), qdims=(1,1,1), qstep=0.1)
	
    Nh = qdims[1]; Nk = qdims[2]; Nl = qdims[3]; Δq = qstep
    Qs = [[h,k,l] for h in -Nh:Δq:Nh, k in -Nk:Δq:Nk, l in -Nl:Δq:Nl]
	
    swt2D = SpinWaveTheory(sys; measure=ssf_perp(sys))
    result2D = intensities(swt2D, Qs[:]; energies, kernel=kern)
    
    #This section is still suspect. 
    #Need a good path sampling scheme.
    #Might move this section to it's own function.
    @inline function fibonacci_sphere(n::Int)
        ϕ = (1 + sqrt(5))/2; ℓ = 0:(n-1)
        zs = @. 1 - 2ℓ/(n-1); rs = @. sqrt(1 - zs^2); θs = @. 2π*ℓ/ϕ; 
        xs = @. rs * cos(θs); ys = @. rs .* sin(θs); 
        [[xs[i], ys[i], zs[i]]/sqrt(xs[i]^2 + ys[i]^2 + zs[i]^2) for i in 1:n]
    end
    @inline function qline(u::AbstractVector{<:Real}; qmax=1.0, nhkl=701)
        t = range(-qmax, qmax; length=nhkl)
        [[t[i]*u[1], t[i]*u[2], t[i]*u[3]] for i in eachindex(t)]
    end
    function generate_paths(cryst; qmax=1.0, nhkl=701, n_dir=96, 
                                    include_axes=true, 
                                    include_face_diagonals=true,
                                    include_body_diagonals=true)
        dirs = Vector{Vector{Float64}}()
        if include_axes
            append!(dirs, [[1,0,0], [0,1,0], [0,0,1], [-1,0,0], [0,-1,0], [0,0,-1]])
        end
        if include_face_diagonals
            append!(dirs, [[1,1,0], [1,0,1], [0,1,1], -[1,1,0], -[1,0,1], -[0,1,1]])
        end
        if include_body_diagonals
            append!(dirs, [[1,1,1], -[1,1,1]])
        end
        for u in fibonacci_sphere(n_dir)
            push!(dirs, collect(u))
        end
        #Normalize & Remove Degeneracy
        dirs = unique!(map(u -> (n = sqrt(sum(abs2, u)); [u[1]/n, u[2]/n, u[3]/n]), dirs))
        return [q_space_path(cryst, qline(u; qmax=qmax, nhkl=nhkl), nhkl) for u in dirs]
    end
    paths = generate_paths(cryst; qmax=1.0, nhkl=101, n_dir=36, 
                            include_axes=true, include_face_diagonals=true, include_body_diagonals=true)


    #Could probably organize the output better. 
    results = Vector{Sunny.Intensities{Float64, Sunny.QPath, 2}}(undef, length(paths))
    @threads for i in eachindex(paths)
	    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
	    results[i] = intensities(swt, paths[i]; energies, kernel=kern)
    end
    return [result2D, Qs], [results, paths], energies
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
