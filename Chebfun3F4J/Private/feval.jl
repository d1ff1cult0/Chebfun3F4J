using ApproxFun
using LinearAlgebra

function feval(cf3F, x :: Float64, y :: Float64, z :: Float64)
    r1 = cf3F.r[1]
    r2 = cf3F.r[2]
    r3 = cf3F.r[3]
    
    u = interpolateColumns(cf3F.U,x,r1)
    v = interpolateColumns(cf3F.V,y,r2)
    w = interpolateColumns(cf3F.W,z,r3)

    out = cf3F.C
    out = u'*reshape(out, r1, r2*r3)
    out = v'*reshape(out, r2, r3)*w
    
    return out
end

function interpolateColumns(U, x, r)
    res = Vector{Float64}(undef, r)
    for i in eachindex(U)
        res[i] = U[i](x)
    end
    return res
end