module SLHashes

using LinearAlgebra: Tridiagonal, I
export get_slhash


function get_slhash(first_lambda::Int, n::Int, a::Int, b::Int, l::Int, A_map::Vector{Int}, Ainv_map::Vector{Int}, B_map::Vector{Int}, Binv_map::Vector{Int})
	A::Matrix{Int} = Tridiagonal(zeros(n-1), ones(n), a*ones(n-1))^l
	B::Matrix{Int} = Tridiagonal(b*ones(n-1), ones(n), zeros(n-1))^l
	Ainv::Matrix{Int} = A^-1
	Binv::Matrix{Int} = B^-1
	matrices = (A, B, Ainv, Binv)

	mappings = (Tuple(A_map), Tuple(B_map), Tuple(Ainv_map), Tuple(Binv_map))

	function slhash(sequence::Vector)::Matrix{Int}
		lambda = first_lambda
		result = Matrix{Int}(I, n, n)

		for x in sequence
			choice = mappings[lambda][x]
			result *= matrices[choice]
			lambda = (choice + 2) % 4
		end
		return result
	end
end


end
