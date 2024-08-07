module SLHashes

using LinearAlgebra: Tridiagonal, I
using Primes: isprime, factor
using AbstractAlgebra: GF, identity_matrix, matrix

export get_slhash, get_mapping


include("hash.jl")
include("crack.jl")


end
