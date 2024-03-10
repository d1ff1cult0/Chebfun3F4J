include("constructor.jl")

function basictestfun(x,y,z)
    return 1 ./ (1 .+ 25 * (x.^2 + y.^2 + z.^2))
end

cf3F(basictestfun)