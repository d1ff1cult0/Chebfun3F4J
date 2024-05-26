include("Chebfun3F4J/constructor.jl")
include("Chebfun3F4J/phase1.jl")
include("testfunction.jl")

using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 15 #Stabilizeerd de benchmark (kinda ig idrk)
function execute_test()
    for i in 42:42 
        try
            f = testfunction(i)
            cf = cf3F(f)
        catch
            println("fail detected")
        end
    end
end

function testPhase1()
    f = testfunction(42)
    n = 100

    a = 0
    b = 0
    c = 0
    for _ in 1:n
        t1, t2, t3 = phase1(f)
        a += t1
        b += t2
        c += t3
    end
    return a, b, c
end

#cf = cf3F(testfunction(1), true); println("sd")
#@benchmark execute_test()
testPhase1()