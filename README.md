# SLHashes

[![Build Status](https://github.com/ZeeZedZhi/SLHashes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ZeeZedZhi/SLHashes.jl/actions/workflows/CI.yml?query=branch%3Amain)

An implementation of the family of hash functions described in [Post-quantum hash functions using special linear groups over finite fields](https://eprint.iacr.org/2022/896).

## Specifying Hash Function Parameters and Usage of Hash Function
Example of ```get_slhash``` usage. ```n```, ```a```, ```b```, and ```l``` are as described in Definition 2.4 of the paper. ```p``` is an optional parameter; if it is not prime, then it is ignored.
```
# Specifying hash function paremeters
first_lambda = 1
n = 3
a = 4
b = 2
l = 4
p = 0

choices = [2 3 4; 1 3 4; 1 4 2; 1 3 2]

slhash = get_slhash(first_lambda, n, a, b, l, choices, p)

# Hashing the sequence 2232221
hash = slhash([2, 2, 3, 2, 2, 2, 1]) #returns [694190977 233260720 29297952; -38379648 -12896255 -1619792; 1191936 400512 50305]
```

```choices``` is a 4x3 matrix to specify $s_\lambda$. In accordance with the paper's definition, ```get_slhash``` is set up such that $s(1) = A$, $s(2) = B$, $s(3) = A^{-1}$, $s(4) = B^{-1}. Thus, the matrix 
```math
  \begin{pmatrix}
    2 & 3 & 4\\
    1 & 3 & 4\\
    1 & 4 & 2\\
    1 & 3 & 2
  \end{pmatrix}
```
corresponds with the example given in Definition 2.6.
