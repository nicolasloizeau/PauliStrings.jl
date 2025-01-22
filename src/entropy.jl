

function renyi_entropy(o::Operator, alpha::Int)
    o = compress(o)
    c = get_coefs(o)
    c /= sum(abs.(c))
    return log(sum(abs.(c) .^ alpha)) / (1 - alpha)
end
