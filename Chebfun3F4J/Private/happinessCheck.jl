using ApproxFun
function happinessCheck(coeffs, tol)
    f = Fun(Chebyshev(), vec(coeffs))
    fc = ApproxFun.chop(f, tol) 
    return length(fc.coefficients) < length(coeffs)
end