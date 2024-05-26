#TODO: rekening houden met imaginaire functies in chebtech2, zie valls2coeffs in matlab
using FFTW
using FastChebInterp
using BasicInterpolators

function createCT2(W :: Matrix{Float64})
    coeffs = chebtech2(W)
    coeffs = sum(abs.(coeffs), dims=2)
    return coeffs
end

#decommissioned
function chebfunOld(W :: Matrix{Float64}, tol::Float64=eps())
    P = Vector{FastChebInterp.ChebPoly{1, Float64, Int64}}(undef, size(W,2)) # Preallocate memory
    for i in 1:size(W,2)
        col = W[:,i]
        F = chebinterp(col, -1, 1,tol=tol)
        P[i] = F
    end
    return P
end

function chebfun(W :: Matrix{Float64})
    P = Vector{ChebyshevInterpolator{size(W,1), Float64}}(undef, size(W,2)) # Preallocate memory
    grid = chebygrid(size(W,1))
    for i in 1:size(W,2)
        col = W[:,i]
        F = ChebyshevInterpolator(grid,col)
        P[i] = F
    end
    return P
end

function chebtech2(W :: Matrix{Float64})
    n = size(W,1)
    tmp = [W[n:-1:2,:] ; W[1:n-1,:]]
    coeffs = real(ifft(tmp, (1,)))
    coeffs = coeffs[1:n,:]
    coeffs[2:n-1,:] = 2*coeffs[2:n-1,:]
    return coeffs
end