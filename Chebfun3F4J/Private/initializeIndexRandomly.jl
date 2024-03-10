function initializeIndexRandomly(r, maxVal)
    box = floor(maxVal/r)
    X = []
    for i in 1:r
        push!(X, i*box + rand(1:box))
    end
    return X
end

function initializePresetIndex(r, maxVal)
    if r == 6 && maxVal == 17
        J = [4, 6, 8, 9, 12, 14]
        K = [4, 5, 8, 10, 11, 14]
    end
    if r == 9 && maxVal == 33
        J = [4, 8, 11, 14, 16, 19, 22, 27, 29]
        K = [6, 7, 12, 15, 17, 20, 22, 27, 28]
    end
    if r == 14 && maxVal == 46
        J = [5,9,12,14,16,21,22,27,28,32,24,28,40,43]
        K = [6,8,11,14,17,19,24,27,30,33,35,38,41,44]
    end
    return J,K
end