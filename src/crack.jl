function get_mapping(slhash::Function, a::Int, b::Int, l::Int, p::Int=0)::Matrix{Int}
	n = size(slhash([1]), 1)


	A::Matrix{Int} = Tridiagonal(zeros(n-1), ones(n), a*ones(n-1))^l
	B::Matrix{Int} = Tridiagonal(b*ones(n-1), ones(n), zeros(n-1))^l
	Ainv::Matrix{Int} = A^-1
	Binv::Matrix{Int} = B^-1

	ISPRIME = false
	if isprime(p)
		ISPRIME = true
		A = mod.(A, p)
		B = mod.(B, p)
		Ainv = mod.(Ainv, p)
		Binv = mod.(Binv, p)
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
				mappings[inversion[mappings[X, x]], y] = matrixToLambda[ISPRIME ? mod.(matrix, p) : matrix]
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
						mappings[inversion[X], z] = matrixToLambda[ISPRIME ? mod.(matrix, p) : matrix]
					end
					break
				end
			end
			break
		end
	end

	return mappings
end