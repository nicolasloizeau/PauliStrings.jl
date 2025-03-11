
"""
convert a string in the computational basis to a state vector
"""
function state_to_vect(state::String)
    v = 1
    for i in 1:length(state)
        if state[i] == '0'
            v = kron(v, [1, 0])
        else
            v = kron(v, [0, 1])
        end
    end
    return v
end


@testset "states" begin
    Z = Operator(1)
    Z += "Z"
    @test expect(Z, "0") == 1
    @test expect(Z, "1") == -1
    N = 4
    O1 = rand_local2(N) + rand_local1(N)
    O2 = rand_local2(N) + rand_local1(N)
    O1dense = op_to_dense(O1)
    state = join(rand(["0","1"], N))
    statevect = state_to_vect(state)
    @test abs(expect(O1, state) - statevect' * O1dense * statevect)<1e-10
    @test abs(expect_product(O1, O2, state)-expect(O1*O2, state))<1e-10
end
