push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots

L = 10
r = 1
V = 12.4

f(h) = r * sin(0.5π - h/r^2*sqrt(r^2 - h^2) - V/(L*r)) - h
c, itt = smartBisectionMethod(f, 10., 50.; domain = (-1.,1.))

x = LinRange(-1,1,9000)
y = f.(x)
plot(x, y, label = "Function", xlabel = "Height", ylabel = "f(h)", title = "Function f(h) = L*(0.5π*r^2 -r^2*asin(h/r) - h*sqrt(r^2 - h^2)) - V", lw = 2, color = :red)
plot!(x, zeros(length(x)), label = "f(h) = 0", lw = 2, color = :blue)
plot!(ones(length(x)) .* c, LinRange(minimum(y), maximum(y), length(y)), label = "x = $c", lw = 2, color = :blue)
scatter!([c], [0], label = "Root", color = :blue)