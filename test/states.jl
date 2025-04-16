
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
    state1 = join(rand(["0","1"], N))
    statevect1 = state_to_vect(state1)
    state2 = join(rand(["0","1"], N))
    statevect2 = state_to_vect(state2)
    @test abs(expect(O1, state1) - statevect1' * O1dense * statevect1)<1e-10
    @test abs(expect(O1, state1, state2) - statevect2' * O1dense * statevect1)<1e-10
    @test abs(expect_product(O1, O2, state1)-expect(O1*O2, state1))<1e-10

end
