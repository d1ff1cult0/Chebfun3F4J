include("Chebfun3F4J/constructor.jl")
include("Test/testfunction.jl")
include("Chebfun3F4J/private/sampleTest.jl")

using Statistics

function getWorstAndAverageError(testfuncNumber)
    maximums = []
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
        push!(maximums, maximum(map(abs, vHandle-vFun)))
        push!(averages, mean(map(abs, vHandle-vFun)))
    end

    

    println("Average of maximum errors in Halton Points: error $(mean(maximums))")
    println("Average of average errors in Halton Points: error $(mean(averages))")
end
