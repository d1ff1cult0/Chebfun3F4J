# Adaptive Cross Approximation with full pivoting
function ACA(A, tol, maxIter)
    Ac = []; Ar = []; At = []
    rowInd = []; colInd = []
    Aoriginal = deepcopy(A)

    for iter in 1:maxIter
        error, I2 = findmax(map(abs, A))
        if isempty(error) || error < tol
            Ac = Aoriginal[:,colInd]
            At = Aoriginal[rowInd,colInd]
            Ar = Aoriginal[rowInd,:]
            return [Ac, At, Ar, rowInd, colInd]
        end
        I = I2[1]; J = I2[2]
        push!(rowInd, I)
        push!(colInd, J)
        A = A .- (A[:,J] * A[I,:]' ./ A[I,J])
    end
    Ac = Aoriginal[:, colInd]
    Ar = Aoriginal[rowInd, :]
    At = Aoriginal[rowInd, colInd]
    return [Ac, At, Ar, rowInd, colInd]
end