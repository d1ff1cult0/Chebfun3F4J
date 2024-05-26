function testfunction(i)
    if i == 1 # Werkt
        return (x,y,z) -> 1 ./ (1 .+ 25 * (x.^2 + y.^2 + z.^2))
    elseif i == 2 # Werkt
        return (x,y,z) -> exp.(x.*y.*z)
    elseif i == 3 # Werkt niet -> stackoverflow
        return (x,y,z) -> 0
    elseif i == 4 # Werkt
        return (x,y,z) -> 3*(1 .- x).^2 .* exp.(-x.^2 .- (y .+ 1).^2) - 10*(x/5 .- x.^3 .- y.^5) .* exp.(-x.^2 .- y.^2) - (1/3) * exp.(-(x .+ 1).^2 .- y.^2)
    elseif i == 5 # Werkt
        return (x,y,z) -> -cos.(50*pi*(x+y+z))
    elseif i == 6 # Werkt
        return (x,y,z) -> exp.(-(x-y))
    elseif i == 7 # Werkt
        return (x,y,z) -> 1 ./ cosh.(3*(x+y+z)).^2
    elseif i == 8 # Werkt
        return (x,y,z) -> tanh.(5*(x+z)) .* exp.(y)
    elseif i == 9 # Werkt
        return (x,y,z) -> (x.^2 .- y.^3 .+ 1/8)
    elseif i == 10 # Werkt
        return (x,y,z) -> x.^10 .* y.^10 .* z.^10
    elseif i == 11 # Werkt niet, wel in matlab -> Stackoverflow
        return (x,y,z) -> 1
    elseif i == 12 # Werkt
        return (x,y,z) -> -cos.(50*pi*(x+y+z))
    elseif i == 13 # Werkt
        return (x,y,z) -> x .* y .+ cos.(x.*y)
    elseif i == 14 # Blijft hangen in fase 2 (maar ook in matlab dus ok)
        return (x,y,z) ->  1 ./(0.00001.+x.^2)
    elseif i == 15 # Werkt
        return (x,y,z) -> -1000 .+ x .* y .* z.^100
    elseif i == 16 # Werkt
        return (x,y,z) -> y .* z .+ cos.(z.*y)
    elseif i == 17 # Werkt
        t = 5
        return (x,y,z) -> exp.(-((x.-0.6*cos.(pi*t)).^2 .+ (y.-0.6*sin.(pi*t)).^2)./0.1)
    elseif i == 18 # Werkt
        c = (sqrt(5) + 1) / 2
        return (x,y,z) -> 8 .*(x.^2 .- c.^4 .*y.^2) .*(y.^2 .- c.^4 .*z.^2) .*(z.^2 .- c.^4 .*x.^2) .*(x.^4 .+ y.^4 .+ z.^4 .- 2 .*x.^2 .*y.^2 .- 2 .*x.^2 .*z.^2 .- 2 .*y.^2 .*z.^2) .+ (3 .+ 5 .*c) .* ((x.^2 .+ y.^2 .+ z.^2)).^2 .* ((x.^2 .+ y.^2 .+ z.^2 .- (2 .- c))).^2
    elseif i == 19
        # /
    elseif i == 20 # Werkt
        return (x,y,z) -> ((x.-0.4).^2 .+ y.^2) .* ((x.+0.4).^2 .+ y.^2) - z.^4
    elseif i == 21 # Werkt
        return (x,y,z) -> 81 .*(x.^3 .+ y.^3 .+ z.^3) .- 189 .*(x.^2 .*y .+ x.^2 .*z .+ y.^2 .*x .+ y.^2 .*z .+ z.^2 .*x .+ z.^2 .*y) .+ 54 .*(x.*y.*z) .+ 126 .*(x.*y .+ x.*z .+ y.*z) .- 9 .*(x.^2 .+ y.^2 .+ z.^2) .- 9 .*(x .+ y .+ z) .+ 1
    elseif i == 22 # Werkt niet, complexe getallen
        return (x,y,z) -> abs.(((x .+ 1i.*y).*exp.(1i.*z)).^2 .- 1).^2
    elseif i == 23 # Werkt
        return (x,y,z) -> cos.(pi .+ 5 .* (x+y+z))
    elseif i == 24 # Werkt
        return (x,y,z) -> 1 ./ (16 .+ 5 .* (x+y+z))
    elseif i == 25 # Werkt
        return (x,y,z) -> exp.(-25*(x.^2+y.^2+z.^2))
    elseif i == 26 # Werkt
        return (x,y,z) -> 1 ./ ((0.04 .+ x.^2) .* (0.04 .+ y.^2) .* (0.04 .+ z.^2))
    elseif i == 27 # Werkt
        return (x,y,z) -> x.^2 ./4000 .+ y.^2 ./ 4000 .+ z.^2 ./ 4000 .- cos.(x).*cos.(y ./ sqrt.(2)).*cos.(z ./ sqrt.(3)) .+ 1
    elseif i == 28 # Werkt
        return (x,y,z) -> z.^2 .+ 0.5.*(x.^2 .+ y.^2 .- .5).^2
    elseif i == 29 # Werkt
        return (x,y,z) ->  x.^4 + y.^4 + z.^4 .- 0.1.*(x.^2 .+y.^2 .+ z.^2) .- 0.5.*(x.^2 .* y.^2 .+ x.^2 .* z.^2 .+ y.^2 .* z.^2) .+ 0.1*x.*y.*z .- 1
    elseif i == 30 # Werkt
        return (x,y,z) -> cos.(2 .* pi .* x).^2 .+ cos.(2 .*pi.*y).^2 .+ cos.(2 .* pi .*z).^2
    elseif i == 31 # Werkt
        return (x,y,z) -> (sin.(pi.*(1 .+ (x .- 1)./4))).^2 .+ ((z .- 1)./4).^2 .* (1 .+(sin.(2 .*pi.*(z .- 1)./4)).^2) .+ (((x .- 1)./4)).^2 .* (1 .+10 .*(sin.(pi.*(1 .+ (x .- 1)./4).+1)).^2) .+ (((y .- 1)./4)).^2 .* (1 .+10 .*(sin.(pi.*(1 .+ (y .- 1)./4).+1)).^2)
    # pas op, vanaf hier komen de functies niet meer overeen met de functies uit de matlab code aangezien functie 31 2x gedefinieerd is 
    elseif i == 32 # Werkt
        return (x,y,z) -> 100 ./(1 .+100 .*(x.^2+y.^2+z.^2))
    elseif i == 33 # Werkt
        return (x,y,z) -> 30 .+ (x.^2 .- 10 .*cos.(2 .*pi.*x)) .+ (y.^2 .- 10 .*cos.(2 .*pi.*y)) .+ (z.^2 .- 10 .*cos.(2 .*pi.*z))
    elseif i == 34 # Werkt
        return (x,y,z) -> 100 .*(y .- x.^2).^2 .+ (x .- 1).^2 .+ 100 .*(z .- y.^2).^2 .+ (y .- 1).^2
    elseif i == 35 # Werkt
        return (x,y,z) -> 1 ./(1 .+x.^2 .+y.^2 .+z.^2)
    elseif i == 36 # Werkt
        return (x,y,z) -> (cos.(2 .*x.+1) .+ 2 .*cos.(3 .*x.+2) .+ 3 .*cos.(4 .*x.+3) .+ 4 .*cos.(5 .*x.+4) .+ 5 .*cos.(4 .*x.+5)) .* (cos.(2 .*y.+1) .+ 2 .*cos.(3 .*y.+2) .+ 3 .*cos.(4 .*y.+3) .+ 4 .*cos.(5 .*y.+4) .+ 5 .*cos.(4 .*y.+5)) .* (cos.(2 .*z.+1) .+ 2 .*cos.(3 .*z.+2) .+ 3 .*cos.(4 .*z.+3) .+ 4 .*cos.(5 .*z.+4) .+ 5 .*cos.(4 .*z.+5))
    elseif i == 37 # Werkt
        return (x,y,z) -> 3 .*x.^7 .*z .+ y.*z .+ y.*z.^2 .+ log.(2 .+ y).*z.^3 .- 2 .*z.^5
    elseif i == 38 # Werkt
        return (x,y,z) -> exp.(sin.(50 .*x)) .+ sin.(60 .*exp.(y)).*sin.(60 .*z) .+ sin.(70 .*sin.(x)).*cos.(10 .*z) .+ sin.(sin.(80 .*y)) .- sin.(10 .*(x .+ z)) .+ (x.^2 .+ y.^2 .+ z.^2)./4
    elseif i == 39 # Werkt
        return (x,y,z) -> log.(1 .+x.^2 .+y.^2 .+z.^2)
    elseif i == 40 # Werkt
        return (x,y,z) -> y .*z .+ cos.(z.*y)
    elseif i == 41 # Werkt
        return (x,y,z) -> 10000 ./(1 .+10000 .*(x.^2 .+y.^2 .+z.^2))
    elseif i == 42 # Werkt
        return (x,y,z) -> log.(x.+y.*z.+exp.(x.*y.*z).+cos.(sin.(exp.(x.*y.*z))))
    elseif i == 43 # Werkt
        return (x,y,z) -> cos.(2 .*pi.*x).^2 + cos.(2 .*pi .*y).^2 + cos.(2 .*pi .*z).^2
    elseif i == 44 # Werkt
        return (x,y,z) -> exp.(sin.(50 .*x)) .+ sin.(60 .* exp.(y)) .* sin.(60 .*y) .+ sin.(70 .*sin.(x)).*cos.(10 .*z) .+ sin.(sin.(80 .*y)) .- sin.(10 .*(x .+ z)) .+ (x.^2 .+ y.^2 .+ z.^2)./4
    elseif i == 45 # Werkt
        return (x,y,z) -> x.*z .+ x.^2 .*y
    elseif i == 46 # Werkt
        return (x,y,z) -> sin.(1 ./((0.04.+x.^2)).*(0.04.+y.^2).*(0.04.+z.^2))
    elseif i == 47 # Complexe getallen
        #
    elseif i == 48 # Werkt
        return (x,y,z) -> 1 ./(cosh.(3 .*(x.+y.+z)).^2)
    elseif i == 49 # Werkt
        return (x,y,z) -> exp.(sin.(x.+2 .*y.+3 .*z)).+y .*z
    elseif i == 50
        # Complexe getallen
    elseif i == 51
        # Complexe getallen
    elseif i == 52 # Werkt
        return (x,y,z) -> 1 ./((x.+1.1).+(y.+1.1).+(z.+1.1))
    elseif i == 53 # Max number of restarts
        return (x,y,z) -> exp.(-sqrt.((x.-1).^2 .+(y.-1).^2 .+ (z.-1).^2))
    elseif i == 54
        # /
    else
        return (x,y,z) ->  0
    end
end