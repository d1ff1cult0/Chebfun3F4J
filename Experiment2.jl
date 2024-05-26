# Test influence of grid size on speed and accuracy
include("Chebfun3F4J/constructor.jl")
include("Test/exp2Constructor.jl")
include("Test/testfunction.jl")

using Statistics
using Plots
using Unzip

function getSpeedByGridSize(testfuncNumber :: Int64)
    timOfGridSize = []
    for g in 10:10:250
        times = []
        for k in 1:50
            f = testfunction(testfuncNumber)
            t = time()
            cf = cf3F(f,false,false,g)
            push!(times, time() - t)
            cf = nothing
        end
        tim = mean(times)
        println("grid size $(g) done")
        push!(timOfGridSize, (g,tim))
    end 
    z = unzip(timOfGridSize)
    for i in eachindex(z[1])
        println(z[1][i])
    end
    for i in eachindex(z[2])
        println(z[2][i])
    end
    p = plot(z[1],z[2], label="Average time")
    title!("Speed of grid sizes")
    xlabel!("g")
    ylabel!("Time (s)")
    savefig("SpeedByGridSizeFunction$(testfuncNumber).pdf")
end

function getErrorByGridSize(testfuncNumber :: Int64)
    f = testfunction(testfuncNumber)
    cf, data = cf3F2(f)
    z = unzip(data)
    z1 = [minimum(a) for a in z[1]]
    println(z1)
    println(z[2])
    println(z[3])
    p = plot(z1, log10.(z[2]),label = "Max Error")
    p = plot!(p,z1, log10.(z[3]),label = "Average Error")
    title!("Accuracy of grid sizes")
    xlabel!("g")
    ylabel!("Error")
    cf = nothing
    savefig("getErrorByGridSize$(testfuncNumber).pdf")
end

function getExp2Data()
    println("Data collection started")
    for i in 1:53
        try
            if i <= 20 # don't test every function
                getSpeedByGridSize(i)
            end
            if i in [32,41,46] # only functions with restarts >= 2
                getErrorByGridSize(i)
            end
            println("function $(i) done")
        catch
            println("function $(i) failed")
        end
    end
end

getSpeedByGridSize(13)