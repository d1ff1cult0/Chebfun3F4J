function evalTensor(I, J, K, ff, vectorize)
    if vectorize 
        n = [numel(I), numel(J), numel(K)]
        x = zeros(Int64, n[1], 1, 1)
        x[:, 1, 1] = I
        X = repeat(x, 1, n[2], n[3])
        y = zeros(Int64, 1, n[2], 1)
        y[1, :, 1] = J
        Y = repeat(y, n[1], 1, n[3])
        z = zeros(Int64, 1, 1, n[3])
        z[1, 1, :] = K
        Z = repeat(z, n[1], n[2], 1)       
        return ff(X,Y,Z)
    else
        T = zeros(Float64, size(I, 2), size(J,2), size(K, 2))
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

function numel(A)
    return prod(size(A))
end