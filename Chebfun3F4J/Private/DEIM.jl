function DEIM(U::Matrix{Float64})
    indices = []
    I = findmax(map(abs, U[:,1]))[2]
    push!(indices, I)
    for l in 2:size(U,2)
        c = U[indices, 1:(l-1)] \ U[indices, l]
        r = U[:, l] - U[:, 1:(l-1)]*c
        I = findmax(map(abs, r))[2][1]
        push!(indices, I)
    end
    return Int64.(indices)
end