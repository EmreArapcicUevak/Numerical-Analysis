push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

#### Other modules that you might need go here ####

###################################################
using LinearSystemModule
#### Your code goes here ####

A = [1 0 0 0; 1 π/2 (π/2)^2 (π/2)^3; 1 π (π)^2 (π)^3; 1 3π/2 (3π/2)^2 (3π/2)^3]
b = Float64[0, 1, 0, -1]

result = PLUFactorization(A)
coef = SolvePLU(result, b)

x = LinRange(0, 6, 9000)
y = sin.(x)
poly = coef[1] .+ coef[2] .* x .+ coef[3] .* x.^2 .+ coef[4] .* x.^3

plot(x, y, label = "sin(x)")
scatter!([0, pi/2, pi, 3pi/2],[0,1,0,-1], label = "Value Points")
plot!(x, poly, label = "Interpolated Polynomial")