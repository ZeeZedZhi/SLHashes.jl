using SLHashes
using Test

@testset "SLHashes.jl" begin
    
    sequence = [2, 2, 3, 2, 2, 2, 1]
	testhash1 = get_slhash(first_lambda=1, n=3, a=4, b=2, l=4, mappings=[2 3 4; 1 3 4; 1 4 2; 1 3 2])
	@test testhash1(sequence) = [694190977 233260720 29297952; -38379648 -12896255 -1619792; 1191936 400512 50305]

end
