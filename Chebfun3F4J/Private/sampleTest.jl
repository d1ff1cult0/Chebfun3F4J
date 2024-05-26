using HaltonSequences
include("feval.jl")

function sampleTest(cf3F, f :: Function, tol :: Float64, enableLog :: Bool)
    n = 30

    (x, y, z) = halton(n)
    vFun = zeros(Float64, n); vHandle = zeros(Float64, n)
    for i in 1:n
        vFun[i] = f(x[i],y[i],z[i])
        vHandle[i] = feval(cf3F,x[i],y[i],z[i])
    end

    diff = abs.(vHandle .- vFun)
    maxError = maximum(diff)
    meanError = mean(diff)
    tolCheck = maxError <= 10 * tol

    if enableLog
        println("Sample Test: error $(maxError) tol $(10*tol)")
    end

    return tolCheck, maxError, meanError
end

function halton(numpts)
    H = HaltonPoint(3)[1:numpts]
    x = 2 .* getindex.(H, 1) .- 1
    y = 2 .* getindex.(H, 2) .- 1
    z = 2 .* getindex.(H, 3) .- 1
    return x, y, z
end