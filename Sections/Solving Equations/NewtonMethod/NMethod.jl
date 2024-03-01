push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = x * exp(-x)
f_prime(x) = exp(-x) - x * exp(-x)

result = NewtonMethod(f, f_prime, 2.5; maxIterations = 100, ϵ = 1e-12)
result


f(x) = x^5
f_prime(x) = 5x^4

result = NewtonMethod(f, f_prime, 1)
ShowPointModule.showPoint(f,(result.c,0); domain = (-20, 20))

f(x) = sin(x)
f_prime(x) = cos(x)
result1 = NewtonMethod(f, f_prime, 2)

g(x) = sin(x)^2
g_prime(x) = 2sin(x)cos(x)
result2 = NewtonMethod(g, g_prime, 2)

plot(result1.errors ,xlabel = "Itterations", ylabel = "Error" , label = "sin(x)", lw = 2, color = :red)
plot!(result2.errors , label = "sin(x)^2", lw = 2, color = :blue)

f(x) = atan(x)
f_prime(x) = 1 / (1 + x^2)
result3 = NewtonMethod(f, f_prime, 1)
result4 = NewtonMethod(f, f_prime, 2)

x = LinRange(-10, 10, 1000)
plot(x, f.(x))
scatter!([2],[f(2)])

f(x) = 2 - exp(x)
result = QuasiNewtonMethod(f,0,1; ϵ = 1e-12)
ShowPoint(f, result.c; domain = (-2, 2))


g(x) = x^2 - 2
result = QuasiNewtonMethod(g,0,-1; ϵ = 1e-12)
ShowPoint(g, result.c; domain = (-2, 2))

smartBisectionMethod(f, -2, 2)