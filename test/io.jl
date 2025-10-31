


using PauliStrings: paulistringtype


@testset "io" begin
    o1 = construction_example1()
    getproperty.(o1.strings, :w) == [0, 1, 2]
    o2 = construction_example2()
    getproperty.(o2.strings, :v) == [6, 0, 132]
    getproperty.(o2.strings, :w) == [2, 3, 0]
    o3 = construction_example3()
    getproperty.(o3.strings, :v) == [0, 6]
    getproperty.(o3.strings, :w) == [12, 3]
    @test get_coeff(o1, eltype(o1.strings)(4, 0)) == 2
    @test get_coeff(o1, eltype(o1.strings)(0, 1)) == 1
    @test get_coeff(o1, eltype(o1.strings)(2, 2)) == 1
    @test norm(get_coeffs(sort(o1)) - [1, 1, 2]) < 1e-10
    @test norm(get_coeffs(sort(o2)) - [1, 1.5, 2]) < 1e-10
    @test norm(get_coeffs(sort(o3)) - [1, 5]) < 1e-10

    o = Operator(2)
    o += "XY"
    o += 5, "1Z"
    @test op_to_dense(o) == [5.0+0.0im 0.0+0.0im 0.0+0.0im 0.0-1.0im; 0.0+0.0im -5.0+0.0im 0.0+1.0im 0.0+0.0im; 0.0+0.0im 0.0-1.0im 5.0+0.0im 0.0+0.0im; 0.0+1.0im 0.0+0.0im 0.0+0.0im -5.0+0.0im]


    N = 4
    o = rand_local2(N)
    o2 = typeof(o)()
    for i in 1:length(o)
        o2 += Operator(get_pauli(o, i))
    end
    @test opnorm(o2 - all_k_local(N, 2)) < 1e-10
    set_coeffs(o, ones(length(o)))
    @test abs(opnorm(o)^2 - 2^N * length(o)) < 1e-10


    for N in (10, 70)
        for i in 1:10
            O = Operator(N)
            st = randstring(N)
            O += st
            @test string(only(O.strings)) == st
        end
    end
    c, s = op_to_strings(sort(o3))
    @test c == [1.0 - 0.0im, 5.0 + 0.0im]
    @test s == ["XYZ111", "11XX11"]

end



@testset "io single string" begin
    string1 = PauliString("XY111")
    string2 = PauliString{5}(2, 3)
    string3 = PauliString{5}("X", 1, "Y", 2)
    string4 = PauliStringTS{(5,)}("X", 1, "Y", 2)
    string5 = PauliStringTS{(5,)}("XY111")
    @test string1 == string2 == string3
    @test string4 == string5
end
