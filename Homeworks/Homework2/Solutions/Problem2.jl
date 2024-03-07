push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = x^4
f_prime(x) = 4x^3

result = NewtonMethod(f, f_prime, 1)
ShowPoint(f,(result.c,0); domain = (-20, 20))

result2 = QuasiNewtonMethod(f, 0.5, 1)
scatter!([result2.c],[f(result2.c)], label= "Quasi-Newton Method", color = :red, markersize = 5, markerstrokewidth = 0, markerstrokecolor = :red, marker = :circle)
