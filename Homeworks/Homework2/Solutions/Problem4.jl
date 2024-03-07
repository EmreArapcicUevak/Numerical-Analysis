push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

f(x) = sin(x) + x^2 * cos(x) - x^2 - x
df(x) = (1 + 2x) * (cos(x) - 1) -x^2 * sin(x)
d²f(x) = -sin(x) * (1 + 4x) + cos(x) * (2 - x^2) - 2
d³f(x) = sin(x)*(x^2 - 6) - cos(x) * (6x + 1)

f(0)
df(0)
d²f(0)
d³f(0)

# Since the third derivitive is not zero, the multiplicity of the root r = 0 is 3
# This means roots of f, f' and f'' are the same ,but f and f'  converge linearly to r = 0
# But the root of f'' is different and converges quadratically to r = 0

c, itt1 = NewtonMethod(f, df, 1)
c, itt2 = NewtonMethod(df, d²f, 1)
c, itt3 = NewtonMethod(d²f, d³f, 1)

println("Itterations for f: $itt1")
println("Itterations for f': $itt2")
println("Itterations for f'': $itt3")