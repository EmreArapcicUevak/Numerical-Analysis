push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = hornerMethodPolynomial([-3, 1, -2, 1], x)
f_prime(x) = hornerMethodPolynomial([1, -4, 3], x)
result = NewtonMethod(f, f_prime, 2)
ShowPoint(f,(result.c,0); domain = (-20, 20))