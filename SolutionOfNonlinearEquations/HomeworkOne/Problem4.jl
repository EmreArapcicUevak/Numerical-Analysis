push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule

f(x) = abs(x) == 1 ? 5 : x^2 - 1 

result = smartBisectionMethod(f, .5, .6)
c, itt = result
ShowPointModule.showPoint(f, c, domain = (-3, 3))