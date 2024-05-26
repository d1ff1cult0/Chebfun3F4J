
function ind2sub(i1 :: Int64, ndx :: Vector{Int64})
    vi = rem.(ndx .- 1, i1) .+ 1
    v2 = (ndx .- vi) ./ i1 .+ 1
    v1 = vi
    return v1, Int.(v2)
end