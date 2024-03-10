using HaltonSequences
include("feval.jl")

function sampleTest(cf3F, f, tol, enableLog)
    n = 30

    (x, y, z) = halton(n)
    vFun = zeros(n); vHandle = zeros(n)
    for i in 1:n
        vFun[i] = f(x[i],y[i],z[i])
        vHandle[i] = feval(cf3F,x[i],y[i],z[i])
    end

    if enableLog
        println("Sample Test: error $(maximum(map(abs, vHandle-vFun))) tol $(10*tol)")
    end

    return !any(maximum(map(abs, vHandle-vFun)) > 10*tol)
end

function halton(numpts)
    H = HaltonPoint(3)[1:30]
    x = getindex.(H, 1)
    x = 2*x.-1
    y = getindex.(H, 2)
    y = 2*y.-1
    z = getindex.(H, 3)
    z = 2*z.-1
    return [x, y, z]
end