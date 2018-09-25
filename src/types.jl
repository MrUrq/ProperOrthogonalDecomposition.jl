"""
    PODBasis(coefficients::Matrix{Float64}, modes::Matrix{Float64})

Datastructure to store a Proper Orthogonal Decomposition basis.

The original data which the POD basis is representing can be reconstructed by right-
multiplying the modes with the coefficients, i.e. `A = P.modes*P.coefficients`.
"""
struct PODBasis
    coefficients::Matrix{Float64}
    modes::Matrix{Float64}
end
