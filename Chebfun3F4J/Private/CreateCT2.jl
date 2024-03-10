#TODO: rekening houden met imaginaire functies in chebtech2, zie valls2coeffs in matlab (en hoop andere shit die ik overgeslagen heb)
using FFTW

function createCT2(W)
    coeffs = chebtech2(W)
    coeffs = sum(abs.(coeffs), dims=2)
    return coeffs
end

function chebfun(W)
    m = size(W,2)
    vals = zeros(2, m)
    coeffs = chebtech2(W)
    c = deepcopy(coeffs)
    c[2:2:end,:] = -c[2:2:end,:]
    vals[1,:] = sum(c, dims=1)
    vals[2,:] = sum(coeffs, dims=1)
    return vals
end

function chebtech2(W)
    n = size(W,1)
    tmp = [W[n:-1:2,:] ; W[1:n-1,:]]
    coeffs = real(ifft(tmp, (1,)))
    coeffs = coeffs[1:n,:]
    coeffs[2:n-1,:] = 2*coeffs[2:n-1,:]
end