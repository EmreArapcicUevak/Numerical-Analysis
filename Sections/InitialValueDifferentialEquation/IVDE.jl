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

N = 100
f(t :: Real, y :: Real) = -1.2y + 7 * exp(-0.3t)
y₀ = 3  
bounds = (0, 2.5)

solutions = ModifiedEulerMethod(f, y₀, N-1, bounds)
solutions2 = ExplicitEulerMethod(f, y₀, N-1, bounds)
exactSolution(t) = 70/9 * exp(-0.3t) - 43/9 * exp(-1.2t)
points = LinRange(bounds[1], bounds[2], N)

exactPoints = exactSolution.(points)

plot(points, exactPoints, label = "Exact Solution", lw = 2)
plot!(points, solutions, label = "Numerical Solution (Modified)", lw = 2, linestyle = :dash)
plot!(points, solutions2, label = "Numerical Solution (Explicit Euler)", lw = 2, linestyle = :dash)

ErrorModified = log10.(abs.(exactPoints .- solutions))
ErrorExplicit = log10.(abs.(exactPoints .- solutions2))
plot(points, ErrorModified, label = "Error (Modified)", lw = 2)
plot!(points, ErrorExplicit, label = "Error (Explicit Euler)", lw = 2)

boundaries = (0.0, 2*π)
z₀ = Float64[0, 1]
f(t, z) = [z[2], -z[1]]
N = 100

typeof(boundaries)
solutions = ExplicitEulerMethod(f, z₀, N-1, boundaries)
ySolution, yPrimeSolution = Real[], Real[]
for i ∈ solutions
    push!(ySolution, i[1])
    push!(yPrimeSolution, i[2])
end

points = LinRange(boundaries[1], boundaries[2], N)
x(t) = sin(t)

plot(points, x.(points),title = "2nd order differential Equation", label = "Exact Solution", xlabel = "t", ylabel = "y(t)", lw = 3)
plot!(points, ySolution, label = "Y Numerical Solution", lw = 2, linestyle = :dash)
plot!(points, cos.(points), label = "Y' exact solution", lw = 3)
plot!(points, yPrimeSolution, label = "Y' Numerical Solution", lw = 2, linestyle = :dash)

points = LinRange(boundaries[1], boundaries[2], N)




solutions = MidpointMethod(f, z₀, N-1, boundaries)
ySolution, yPrimeSolution = Real[], Real[]
for i ∈ solutions
    push!(ySolution, i[1])
    push!(yPrimeSolution, i[2])
end



points = LinRange(boundaries[1], boundaries[2], N)
plot(points, x.(points),title = "2nd order differential Equation", label = "Exact Solution", xlabel = "t", ylabel = "y(t)", lw = 3)
plot!(points, ySolution, label = "Y Numerical Solution", lw = 2, linestyle = :dash)
plot!(points, cos.(points), label = "Y' exact solution", lw = 3)
plot!(points, yPrimeSolution, label = "Y' Numerical Solution", lw = 2, linestyle = :dash)