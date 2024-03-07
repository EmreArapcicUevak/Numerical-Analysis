push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

F(x :: Vector) = [2(x[1] + 10x[2]) + 40(x[1] - x[4])^3,
    20(x[1] + 10x[2]) + 4(x[3] - 2x[3])^3,
    10(x[3] - x[4]) - 8(x[2] - 2x[3])^3,
    -10(x[3] - x[4]) - 40(x[1] - x[4])^3
]

J(x :: Vector) = [ 
    2 + 120*(x[1] - x[4])^2 20 0 -120*(x[1] - x[4])^2;
    20 200 + 12*(x[2] - 2*x[3])^2 -24*(x[2] - 2*x[3])^2 0;
    0 -24*(x[2] - 2*x[3])^2 10 + 48*(x[2] - 2*x[3])^2 -10;
    -120*(x[1] - x[4])^2 0 -10 10 + 120*(x[1] - x[4])^2
]
c, itt =  MultiDimentionalNewtonMethod(F, J, [3, -1, 0, 1])

F(c)

F(x :: Vector) = [
    exp(x[1]^2 + x[2]^2) - 1,
    exp(x[1]^2 - x[2]^2) - 1
]

J(x :: Vector) = [
    2x[1]*exp(x[1]^2 + x[2]^2) 2x[2]*exp(x[1]^2 + x[2]^2);
    2x[1]*exp(x[1]^2 - x[2]^2) -2x[2]*exp(x[1]^2 - x[2]^2)
]

c, itt =  MultiDimentionalNewtonMethod(F, J, [.1, .1])
F(c)