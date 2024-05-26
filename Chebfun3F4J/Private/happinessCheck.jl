using ApproxFun
function happinessCheck(coeffs :: Matrix{Float64}, tol :: Float64)
    f = Fun(Chebyshev(), vec(coeffs))
    fc = ApproxFun.chop(f, tol) 
    return length(fc.coefficients) < length(coeffs)
end