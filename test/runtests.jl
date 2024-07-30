using SLHashes
using Test

@testset "SLHashes.jl" begin

	n = 3
	a = 4
	b = 2
	l = 4
	mappings = [2 3 4; 1 3 4; 1 4 2; 1 3 2]
    
    testhash = get_slhash(1, n, a, b, l, mappings)
	@test get_mapping(testhash, a, b, l) == mappings

    sequence = [2, 2, 3, 2, 2, 2, 1]
	@test testhash(sequence) == [694190977 233260720 29297952; -38379648 -12896255 -1619792; 1191936 400512 50305]

	primes = [17, 101, 257]
	for prime in primes
		testhashp = get_slhash(1, n, a, b, l, mappings, prime)
		@test get_mapping(testhashp, a, b, l, prime) == mappings
		@test testhashp(sequence) == mod.([694190977 233260720 29297952; -38379648 -12896255 -1619792; 1191936 400512 50305], prime)
	end

end
