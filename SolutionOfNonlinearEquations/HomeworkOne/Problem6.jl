push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = 10x * (1 - ℯ^(-(2e-6)/x)) - 1e-5
x = LinRange(1e-8,1e-4, 9000)
y = f.(x)

plot(x, y, label = "Function", xlabel = "x", ylabel = "f(x)", title = "Function f(x) = 10x * (1 - ℯ^(-2e-6/x)) - 1e-5", lw = 2, color = :red)
plot(x, y .* 1e5, label = "Magnified Function", lw = 2, color = :blue)

result = BisectionMethodModule.smartBisectionMethod(f, 1e-8, 1e-4; domain = (1e-8, 1e-4), precision = 1e-18)
c, itt = result
ShowPointModule.showPoint(f, c, domain = (1e-8, 1e-4))

println("The wanted capacitence is $c F, found in $itt itterations.")