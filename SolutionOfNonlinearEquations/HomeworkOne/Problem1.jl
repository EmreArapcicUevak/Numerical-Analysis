push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots

L = 10
r = 1
V = 12.4

f(h) = L*(0.5π*r^2 -r^2*asin(h/r) - h*sqrt(r^2 - h^2)) - V
println(smartBisectionMethod(f,10.,20.; shift = (x,y) -> .5)) # print's out 1.0
# two roots at -1.0 and 1.0 ignore -1.0 because height can't be negative

x = LinRange(-2,2,9000)
plot(x, f.(x), label = "Function", xlabel = "Height", ylabel = "f(h)", title = "Function f(h) = L*(0.5π*r^2 -r^2*asin(h/r) - h*sqrt(r^2 - h^2)) - V", lw = 2, color = :red)
plot!(x, zeros(length(x)), label = "f(h) = 0", lw = 2, color = :blue)