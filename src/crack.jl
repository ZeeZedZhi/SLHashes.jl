function get_mapping(slhash::Function, a::Int, b::Int, l::Int, p::Int=0)::Matrix{Int}
	n = size(slhash([1]), 1)

	A = Matrix(Tridiagonal(zeros(BigInt, n-1), ones(BigInt, n), fill(BigInt(a), n-1))^l)
	B = Matrix(Tridiagonal(fill(BigInt(b), n-1), ones(BigInt, n), zeros(BigInt, n-1))^l)

	if isprime(p)
		F = GF(p)

		A = matrix(F, A)
		B = matrix(F, B)
		Ainv = A^-1
		Binv = B^-1
	else
		Ainv = BigInt.(A^-1)
		Binv = BigInt.(B^-1)
	end

	matrixToLambda = Dict(A=>1, B=>2, Ainv=>3, Binv=>4)
	mappings::Matrix{Int} = zeros(Int, 4, 3)
	matrices = (A, B, Ainv, Binv)

	# Step 1
	seen = Set(1:4)
	submapping = Dict()
	for j in 1:3
		submapping[j] = matrixToLambda[slhash([j])]
		pop!(seen, submapping[j])
	end

	X = first(seen)
	for j in 1:3
		mappings[X, j] = submapping[j]
	end

	# Step 2
	for x in 1:3
		if mappings[X, x] != inversion[X]
			for y in 1:3
				matrix = matrices[inversion[mappings[X, x]]] * slhash([x, y])
				mappings[inversion[mappings[X, x]], y] = matrixToLambda[matrix]
			end
		end
	end

	# Step 3
	for x in 1:3
		if mappings[X, x] != inversion[X]
			Y = inversion[mappings[X, x]]
			for y in 1:3
				if mappings[Y, y] == X
					for z in 1:3
						matrix = matrices[inversion[X]] * matrices[Y] * slhash([x, y, z])
						mappings[inversion[X], z] = matrixToLambda[matrix]
					end
					break
				end
			end
			break
		end
	end

	return mappings
end