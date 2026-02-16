#CoreUtils.jl (NOT a submodule)
"""Utility macros and helper functions"""

"""
    @name [var1, var2, ...]

Output variable names and values into a named tuple.
Returns `(paths=[var1,var2,...], names=["var1","var2",...])`.

**Note:** Variables must be provided as a literal array expression at the call site.
This macro cannot work with pre-existing arrays.

# Example
```julia
# Correct
H = Sunny.q_space_path(cryst, [[h,0,0] for h in -1:1], 500)
K = Sunny.q_space_path(cryst, [[0,k,0] for k in -1:1], 500)
L = Sunny.q_space_path(cryst, [[0,0,ℓ] for ℓ in -1:1], 500)

qpaths = @name [H,K,L]

#qpaths.paths = [H,K,L]
#qpaths.names = ["H","K","L"]

# Incorrect
arr = [H,K]
qpaths = @name arr #Throws error here
```

"""
macro name(array)
    if !isa(array, Expr) || array.head != :vect
        error("@name requires a literal array expression: @name [H,K,L]")
    end
    elements = array.args
    names = [string(element) for element in elements]
    quote
        paths = [$(esc.(elements)...)]
        names = $names
        (paths=paths,names=names)
    end
end
