push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule

f(x) = 1/x

c, itt = smartBisectionMethod(f, -1 , 10)
ShowPointModule.showPoint(f, c, domain = (-1, 1))