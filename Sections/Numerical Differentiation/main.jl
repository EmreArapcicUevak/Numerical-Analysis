push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

using  NumericalDifferentiationModule

f(x) = sin(x)
fPrime(x) = cos(x)
fPrimePrime(x) = -sin(x)

fAprox = ForwardDifference(f, 1; h = 1e-6)
cos(1)
cos(1) - fAprox

fAprox = CentralDifference(f, 1)
cos(1)
cos(1) - fAprox

fAprox = CentralDifferenceSecondOrder(f, 1)
fPrimePrime(1)
fPrimePrime(1) - fAprox