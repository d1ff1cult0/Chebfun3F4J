function evalTensor(I :: Vector{Int64}, J :: Vector{Int64}, K :: Vector{Int64}, ff, vectorize :: Bool)
    if vectorize 
        return ff(
            repeat(reshape(I, (length(I), 1, 1)), 1, length(J), length(K)),
            repeat(reshape(J, (1, length(J), 1)), length(I), 1, length(K)),
            repeat(reshape(K, (1, 1, length(K))), length(I), length(J), 1)
        )
    else
        T = zeros(Float64, size(I, 2), size(J,2), size(K, 2)) :: Array{Int64,3}
        for i in 1:size(I,2)
            for j in 1:size(J,2)
                for k in 1:size(K,2)
                    T[i, j, k] = ff(I[i], J[j], K[k])
                end
            end
        end
        return T
    end
end
