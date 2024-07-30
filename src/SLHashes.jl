module SLHashes

using LinearAlgebra: Tridiagonal, I
using Primes: isprime, factor

export get_slhash, get_mapping


include("hash.jl")
include("crack.jl")


end
