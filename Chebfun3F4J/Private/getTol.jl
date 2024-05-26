function getTol_old(M :: Matrix{Float64}, pseudoLevel :: Float64, tolOld :: Float64)
    relTol = 2 * size(M,1)^(4/5) * pseudoLevel
    vscale = maximum(abs.(M[:]))
    cheb(i,n) = -cos.((i.-1).*pi/(n-1))
    points = 1:size(M,1)
    points = cheb(points, size(M,1))
    s = 0
    if size(M,1) < size(M,2)
        s = size(M,1)
    else
        s = size(M,2)
    end
    gradNorms = zeros(1,s)
    for i in 1:s
        gradNorms[i] = maximum(abs,(diff(M[:,i]) ./ (diff(points))));
    end
    gradNorms = maximum(gradNorms);
    domDiff = 2;
    absTol = max(gradNorms*domDiff, vscale) * relTol;
    absTol = maximum([absTol, tolOld, pseudoLevel]);
    return [relTol, absTol]
end 

function getTol(M::Matrix{Float64}, pseudoLevel::Float64, tolOld::Float64)
    n = size(M, 1)
    relTol = 2 * n^(4/5) * pseudoLevel
    vscale = LinearAlgebra.norm(M, Inf)
    points = cos.((0:n-1) .* π / (n-1))
    gradNorms = max_abs_diff(M,points)
    domDiff = 2
    absTol = max(gradNorms * domDiff, vscale) * relTol
    absTol = max(absTol, tolOld, pseudoLevel)
    return [relTol, absTol]
end

function max_abs_diff(M, points)
    max_abs_diff_val = zero(eltype(M))
    for i in axes(M,2)
        abs_diff = abs.(diff(M[:, i]) ./ diff(points))
        max_abs_diff_val = max(max_abs_diff_val, maximum(abs_diff))
    end
    return max_abs_diff_val
end