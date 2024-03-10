function feval(cf3F, x, y, z)
    u = cf3F.U[x]
    v = cf3F.V[y]
    w = cf3F.W[z]

    r = rank(cf3F)
    out = cf3F.C
    out = u*reshape(out, r[1], r[2], r[3])
    out = v*reshape(out, r[2], r[3])*w'
    return out
end