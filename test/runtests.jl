using PauliStrings
using Test

N = 10
O = PauliStrings.Operator(N)
O += "X",1


@testset "PauliStrings.jl" begin
    @test PauliStrings.trace(O) == 0
end
