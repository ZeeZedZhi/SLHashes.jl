module SLHashes

using LinearAlgebra: Tridiagonal
export SLMatrices

struct SLMatrices
	A::Matrix{Int}
	B::Matrix{Int}

	function SLMatrices(n::Int, a::Int, b::Int, l::Int)
		return new(Tridiagonal(zeros(n-1), ones(n), a*ones(n-1))^l, Tridiagonal(b*ones(n-1), ones(n), zeros(n-1))^l)
	end
end

end
