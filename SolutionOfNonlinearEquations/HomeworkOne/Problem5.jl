push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = hornerMethodPolynomial([1,-1,-2,1], x)
result = smartBisectionMethod(f, -2, 2)
c, itt = result
ShowPointModule.showPoint(f, c, domain = (-8, 8))

x = LinRange(-8,8,9000)
y₁ = [hornerMethodPolynomial([1,0,-2,1], i) for i ∈ x]
y₂ = x.^2

plot(x, y₁, label = "y₁", xlabel = "x", ylabel = "y(x)", title = "Function f(x) = x^3 - 2x + 1", lw = 2, color = :red)
plot!(x, y₂, label = "y₂", lw = 2, color = :blue)
scatter!([c], [c^3 - 2c + 1], label = "Same Point", color = :blue)