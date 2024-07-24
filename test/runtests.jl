using SLHashes
using Test

@testset "SLHashes.jl" begin
    
    sequence = [2, 2, 3, 2, 2, 2, 1]
	testhash = get_slhash(1, 3, 4, 2, 4, [2 3 4; 1 3 4; 1 4 2; 1 3 2])
	@test testhash(sequence) == [694190977 233260720 29297952; -38379648 -12896255 -1619792; 1191936 400512 50305]


	primes = [17, 101, 257]
	for prime in primes
		testhashp = get_slhash(1, 3, 4, 2, 4, [2 3 4; 1 3 4; 1 4 2; 1 3 2], prime)
		@test testhashp(sequence) == mod.([694190977 233260720 29297952; -38379648 -12896255 -1619792; 1191936 400512 50305], prime)
	end

end
