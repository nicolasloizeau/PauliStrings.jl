


using PauliStrings: paulistringtype

function randstring(N)
    return join(rand(["1", "X", "Y", "Z"], N))
end

@testset "io" begin
    o1 = construction_example1()
    getproperty.(o1.strings, :w) == [0, 1, 2]
    o2 = construction_example2()
    getproperty.(o2.strings, :v) == [6, 0, 132]
    getproperty.(o2.strings, :w) == [2, 3, 0]
    o3 = construction_example3()
    getproperty.(o3.strings, :v) == [0, 6]
    getproperty.(o3.strings, :w) == [12, 3]
    @test get_coef(o1, eltype(o1.strings)(4, 0)) == 2
    @test get_coef(o1, eltype(o1.strings)(0, 1)) == 1
    @test get_coef(o1, eltype(o1.strings)(2, 2)) == 1
    @test norm(sort(get_coeffs(o1)) -[1, 1, 2]) < 1e-10
    @test norm(sort(get_coeffs(o2)) -[1, 1.5, 2]) < 1e-10
    @test norm(sort(get_coeffs(o3)) -[1, 5]) < 1e-10

    println(get_pauli(o1,1))

    for N in (10, 70)
        for i in 1:10
            O = Operator(N)
            st = randstring(N)
            O += st
            @test string(only(O.strings)) == st
            # @test vw_to_string(O.v[1], O.w[1], N)[1] == st
        end
    end




end
