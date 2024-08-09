inversion = (3, 4, 1, 2)


"""
    get_slhash(first_lambda::Int, n::Int, a::Int, b::Int, l::Int, mappings::Matrix{Int}, p::Int=0)

Returns a hash function that maps sequences with terms in {1, 2, 3} to `n` by `n` matrices over the finite field of order `p`. nonprime `p` result in matrices over all integers. 

Use integers to reference matrices; specifically:\n
`1` for `A`\n
`2` for `B`\n
`3` for `A^-1`\n
`4` for `B^-1`\n

Definitions of `n`, `a`, `b`, and `l` are given in the paper. ``first_lambda`` is the matrix that corresponds to the function that maps the first term.

`mappings` is used to choose term mappings as a matrix whose rows respectively correspond to `A`, `B`, `A^-1`, and `B-1` and whose columns respectively correspond to the sequence terms `1`, `2`, and `3`.
"""
function get_slhash(first_lambda::Int, n::Int, a::Int, b::Int, l::Int, mappings::Matrix{Int}, p::Int=0)
	check_parameters(n, a, b, l)
	check_choices(mappings)

	A = Matrix(Tridiagonal(zeros(BigInt, n-1), ones(BigInt, n), fill(BigInt(a), n-1))^l)
	B = Matrix(Tridiagonal(fill(BigInt(b), n-1), ones(BigInt, n), zeros(BigInt, n-1))^l)

	if isprime(p)
		F = GF(p)

		A = matrix(F, A)
		B = matrix(F, B)
		Ainv = A^-1
		Binv = B^-1

		blank = identity_matrix(F, n)
	else
		Ainv = BigInt.(A^-1)
		Binv = BigInt.(B^-1)

		blank = Matrix{BigInt}(I, n, n)
	end

	matrices = (A, B, Ainv, Binv)

	function slhash(sequence::Int...)
		lambda = first_lambda
		result = blank

		for x in sequence
			choice = mappings[lambda, x]
			result *= matrices[choice]
			lambda = inversion[choice]
		end
		return result
	end

	return slhash
end


function check_parameters(n::Int, a::Int, b::Int, l::Int)
	if n == 3
		if mod(a, 3) != 1
			throw(ArgumentError("a is not 1 mod 3"))
		end

		if mod(b, 3) != 2
			throw(ArgumentError("b is not -1 mod 3"))
		end

		if l <= 1 || !isinteger(log(4, l))
			throw(ArgumentError("there is no positive integer k such that 4^k = l"))
		end
	elseif n > 3
		if l < 3*(n-1)
			throw(ArgumentError("l is less than 3(n-1)"))
		end

		if length((local factors = factor(Set, l-1))) != 1 || mod(a, (local q = first(factors))) != 1 || mod(b, q) != 1 || mod(n, q) != 1
			throw(ArgumentError("there is no prime q such that n = a = b = 1 mod q and q^k + 1 = l for some positive integer k"))
		end
	else
		throw(ArgumentError("n is less than 3"))
	end
end


lambdas = Set([1, 2, 3, 4])

function check_choices(mappings::Matrix{Int})
	for i in 1:4
		if (local row = Set(mappings[i, :])) != (local codomain = setdiff(lambdas, i))
			throw(ArgumentError("row $i elements $row does not represent a bijection to $codomain"))
		end
	end
end