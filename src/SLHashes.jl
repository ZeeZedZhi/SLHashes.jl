module SLHashes

using LinearAlgebra: Tridiagonal, I
export get_slhash


function get_slhash(first_lambda::Int, n::Int, a::Int, b::Int, l::Int, mappings::Matrix{Int}, p::Int=0)
	A::Matrix{Int} = Tridiagonal(zeros(n-1), ones(n), a*ones(n-1))^l
	B::Matrix{Int} = Tridiagonal(b*ones(n-1), ones(n), zeros(n-1))^l
	Ainv::Matrix{Int} = A^-1
	Binv::Matrix{Int} = B^-1
	matrices = (A, B, Ainv, Binv)

	if p > 0
		function slhashp(sequence::Vector)::Matrix{Int}
			lambda = first_lambda
			result = Matrix{Int}(I, n, n)

			for x in sequence
				choice = mappings[lambda, x]
				result *= matrices[choice] .% p
				lambda = (choice + 2) % 4
			end
			return result
		end

		return slhashp
	else
		function slhash(sequence::Vector)::Matrix{Int}
			lambda = first_lambda
			result = Matrix{Int}(I, n, n)

			for x in sequence
				choice = mappings[lambda, x]
				result *= matrices[choice]
				lambda = (choice + 2) % 4
			end
			return result
		end

		return slhash
	end
end


end
