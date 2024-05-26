include("Chebfun3F4J/constructor.jl")
include("Test/testfunction.jl")

using Statistics
using Plots

function testFunction(testfuncNumber :: Int64)
    averages = []
    for k in 1:50
        f = testfunction(testfuncNumber)
        cf = cf3F(f)
        
        (x,y,z) = halton(30)
        vFun = zeros(30); vHandle = zeros(30)
        for i in 1:30
            vFun[i] = f(x[i],y[i],z[i])
            vHandle[i] = feval(cf,x[i],y[i],z[i])
        end
        push!(averages, mean(map(abs, vHandle-vFun)))

    end
    return mean(averages)
end

function collectMeans()
    averages = Vector{Float64}(undef,54)
    for i in 1:54
        println(i)
        try
            averages[i] = testFunction(i)
        catch
            averages[i] = 0
        end
    end
    plot(1:54,log10.(averages))
    return averages
end

avgs = collectMeans()
for i in 1:54
    println("$(i): avg = $(avgs[i]) order = $(log10(avgs[i]))")
end
