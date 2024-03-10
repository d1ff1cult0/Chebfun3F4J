
function getTol(M, pseudoLevel, tolOld)
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