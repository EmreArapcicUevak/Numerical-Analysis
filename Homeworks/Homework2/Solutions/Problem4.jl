push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = sin(x) + x^2 * cos(x) - x^2 - x
f_prime(x) = (1 + 2x) * (cos(x) - 1) -x^2 * sin(x)

result = NewtonMethod(f, f_prime, 1)
ShowPoint(f,(result.c,0); domain = (-20, 20))