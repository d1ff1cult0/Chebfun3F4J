include("Chebfun3F4J/Private/sampleTest.jl")
include("Chebfun3F4J/Private/feval.jl")
include("Chebfun3F4J/constructor.jl")
include("Test/testfunction.jl")

using Plots



function experiment4(n :: Int64)
    (x, y, z) = halton(n)

    for functionNumber in 1:n
        try
            f = testfunction(functionNumber)
            

            noiseRange = 10 .^(range(-16,stop=-10,length=70))
            averageErrors = []

            
            averageErrorNoNoise = 0
            
            cfNoNoise = cf3F(f)
            vFun = zeros(n); vHandle = zeros(n)
            for i in 1:n
                vFun[i] = f(x[i],y[i],z[i])
                vHandle[i] = feval(cfNoNoise,x[i],y[i],z[i])
            end
            averageErrorNoNoise = mean(map(abs, vHandle-vFun))
        

            j = 1
            for noiseLevel in noiseRange
                println(j)
                j += 1
                fWithNoise = (x,y,z) -> f(x,y,z) .+ randn()*noiseLevel
                try
                    cf = cf3F(fWithNoise)
                catch
                    @goto escape_label
                end
            

                vFun = []; vHandle = []
                for i in 1:50
                    try
                        push!(vFun, f(x[i],y[i],z[i]))
                        push!(vHandle, feval(cf,x[i],y[i],z[i]))
                    catch
                        println("oopsie i crashed ", j)
                        @goto escape_label
                    end
                end
                push!(averageErrors, mean(abs.(vFun - vHandle)))
            end
            @label escape_label

            p = plot(log10.(noiseRange[1:length(averageErrors)]), log10.(averageErrors),label = "Gemiddelde fout in Halton punten met ruis")
            p = plot!(p, log10.(noiseRange[1:length(averageErrors)]), log10.(fill(averageErrorNoNoise, length(averageErrors))),label = "Gemiddelde fout in Halton punten zonder ruis")
            title!("Evolutie fout in functie van ruis")
            xlabel!("Ruisniveau")
            ylabel!("Gemiddelde fout")
            cf = nothing
            savefig("ExpData4/experiment4-function$(functionNumber).pdf")    
        catch
        end
    end
end