# Adaptive Cross Approximation with full pivoting
function ACA_old(A :: Matrix{Float64}, tol :: Float64, maxIter :: Int64)
    Ac = Vector{Float64}(); Ar = Vector{Float64}(); At = Vector{Float64}()
    rowInd = Vector{Int64}(); colInd = Vector{Int64}()
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
        @views @. A = A - (A[:,J] * A[I,:]' / A[I,J])
    end
    Ac = Aoriginal[:, colInd]
    Ar = Aoriginal[rowInd, :]
    At = Aoriginal[rowInd, colInd]
    return [Ac, At, Ar, rowInd, colInd]
end

function ACA(A::Matrix{Float64}, tol::Float64, maxIter::Int64)
    rowArr = Vector{Int64}(undef, size(A,1))
    colArr = Vector{Int64}(undef, size(A,2))
    Aoriginal = deepcopy(A) :: Matrix{Float64}

    i = 0 :: Int64
    while i <= maxIter 
        i += 1
        error :: Float64, I :: Int64, J :: Int64 = argmax_abs(A)
        if isempty(error) || error < tol
            rowInd = rowArr[1:i-1]
            colInd = colArr[1:i-1]
            return Aoriginal[:,colInd], Aoriginal[rowInd,colInd], Aoriginal[rowInd,:], rowInd, colInd
        end
        rowArr[i], colArr[i] = I, J
        @views @. A .-= A[:, J] * A[I, :]' / A[I, J]
    end
    rowInd = rowArr[1:i]
    colInd = colArr[1:i]
    return [Aoriginal[:,colInd], Aoriginal[rowInd,colInd], Aoriginal[rowInd,:], rowInd, colInd]
end

function argmax_abs(A::Matrix{Float64})
    maxVal = Float64(-1.0)
    maxI :: Int64 = 0
    maxJ :: Int64 = 0
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            absVal = abs(A[i, j]) :: Float64
            if absVal > maxVal
                maxVal = absVal
                maxI, maxJ = i, j
            end
        end
    end
    return maxVal, maxI, maxJ
end