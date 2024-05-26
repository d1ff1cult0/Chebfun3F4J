include("Chebfun3F4J/constructor.jl")
include("Test/testfunction.jl")

using Plots
using Measures

function getAccuracyAtEdgeSolidZ(testfuncNumber)
    f = testfunction(testfuncNumber)
    As = []
    x = Vector(-1:0.1:1)
    y = Vector(-1:0.1:1)
    z = 0.5
    for k in 1:50
        cf = cf3F(f)

        A = Matrix{Float64}(undef,size(x,1),size(y,1))
        A1 = Matrix{Float64}(undef,size(x,1),size(y,1))
        A2 = Matrix{Float64}(undef,size(x,1),size(y,1))
        
        for i in 1:size(x,1)
            for j in 1:size(y,1)
                A1[i,j] = f(x[i],y[j],z)
                A2[i,j] = feval(cf,x[i],y[j],z)
                A = abs.(A1 .- A2)
            end
        end
        push!(As, A)
    end
    # Calculate average error across all 50 samples
    Avgs = Matrix{Float64}(undef,size(x,1),size(y,1))
    for i in 1:size(x,1)
        for j in 1:size(y,1)
            Avg = 0
            for a in As
                Avg += a[i,j]
            end
            Avgs[i,j] = Avg/size(As,1)
        end
    end
    gr(size=(1200,900), dpi=600)
    heatmap(x, y, Avgs,
            colorbar_title="Absolute fout", # You can set fontsize here for colorbar title
            legendfontsize=17, # Adjust legend font size
            titlefontsize=40, # Adjust title font size
            ytickfontsize=17, # Adjust y-axis tick font size
            xtickfontsize=17, # Adjust x-axis tick font size
            yguidefontsize=17, # Adjust y-axis label font size
            xguidefontsize=17, # Adjust x-axis label font size
            legend=:left
    )
    title!("Accuracy of approximation with z=$(z) of function $(testfuncNumber)", fontsize=18) # Adjust title font size
    xlabel!("x", fontsize=20) # Adjust x-axis label font size
    ylabel!("y", fontsize=20) # Adjust y-axis label font size
    savefig("Chebfun3F4J/ExpData3/buh/AccuracyAtEdgeFunction$(testfuncNumber)z$(z).pdf")
end

function getExp3Data(arr)
    for i in arr
        getAccuracyAtEdgeSolidZ(i)

        #try
         #   println("function $(i) done")
        #catch
         #   println("function $(i) fails")
        #end
    end
end