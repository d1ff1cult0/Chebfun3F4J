function initializeIndexRandomly(r :: Int64, maxVal :: Int64)
    box = Int64(floor(maxVal/r))
    X = Vector{Int64}(undef, r)
    for i in 1:r
        X[i] = i*box + rand(1:box)
    end
    return X
end

function initializeIndexPreset(r :: Int64, maxVal :: Int64)
    J = []
    K = []
    if (r == 6 && maxVal == 17)
        J = [4,6,8,9,11,14]
        K = [3,5,7,8,12,14]
    else (r == 12 && maxVal == 46)
        J = [6,7,11,14,18,20,24,26,29,33,34,37]
        K = [4,8,12,15,16,20,23,26,29,32,34,38]
    end
    return J,K
end