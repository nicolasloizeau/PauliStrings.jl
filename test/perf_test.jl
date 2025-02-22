
using PauliStrings
import PauliStrings as ps


function XXZh(N::Int, n, m, h)
    Δ = cos(pi * n / m)
    H = ps.Operator(N)
    H += 'X', 1, 'X', 2
    H += 'Y', 1, 'Y', 2
    H += Δ, 'Z', 1, 'Z', 2
    H += h, 'Z', 1
    return OperatorTS1D(H, full=false)
end


function XXZm3ds1(N::Int)
    H = ps.Operator(N)
    H += -sqrt(3) / 2, "S+", 1, "S+", 2, "S+", 3
    return OperatorTS1D(H, full=false)
end


function get_vmpeak()
    # Only works on Linux systems
    if !isfile("/proc/self/status")
        error("This function only works on Linux systems")
    end
    # Read /proc/self/status and find VmPeak
    for line in eachline("/proc/self/status")
        if startswith(line, "VmPeak:")
            # Extract the number and convert to KB
            return parse(Int, split(strip(line))[2])
        end
    end
    return nothing
end
XXZh2m3(N::Int) = XXZh(N, 2, 3, 2)


N = 40
H = XXZh2m3(N)
O = XXZm3ds1(N)


@time ps.lanczos(H, O, 30, 2^16) # ~ 14.7s
println(get_vmpeak() / 1e6) # ~ 2.6 GB
