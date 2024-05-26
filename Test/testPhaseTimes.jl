include("Chebfun3F4J/constructor.jl")

f = (x,y,z) ->  tanh.(5*(x+z)) .* exp.(y)

t1_times = []
t2_times = []
t3_times = []

for i = 1:100
    handle, t1, t2, t3 = cf3F(f)
    push!(t1_times, t1)
    push!(t2_times, t2)
    push!(t3_times, t3)
end

println("Average t1: ", mean(t1_times))
println("Average t2: ", mean(t2_times))
println("Average t3: ", mean(t3_times))