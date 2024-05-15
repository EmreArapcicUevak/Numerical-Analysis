push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

using DifferentialEquationsModule


f(t, y) = 2y*(1 - y)
y₀ = 1/5
bounds = (10,11)
solutions = ExplicitEulerMethod(f, y₀, 499, bounds)


points = LinRange(bounds[1], bounds[2], 500)
y(t) = 1 / (1 + 4 * exp(2 * (10 - t)))
plot(points, y.(points), label = "Exact Solution", lw = 2)
scatter!(points, solutions, label = "Numerical Solution", lw = 2)