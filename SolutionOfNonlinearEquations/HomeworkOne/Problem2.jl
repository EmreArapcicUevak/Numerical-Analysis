push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots

g = 32.17
f(ω) = 3.4*ω^2 + g*((ℯ^ω - ℯ^(-ω))*0.5 - sin(ω))

c, itt = smartBisectionMethod(f, -5, -0.1; domain = (-Inf, 1e-15), precision = 1e-5)

x = LinRange(-10,10,9000)
y = f.(x)

horLine = zeros(length(x))
verLine = ones(length(x)) .* c

plot(x, y, label = "Function", xlabel = "ω", ylabel = "f(ω)", title = "Function f(ω) = 3.4*ω^2 + g*((ℯ^ω - ℯ^(-ω))*0.5 - sin(ω))", lw = 2, color = :red)
plot!(x, horLine, label = "f(ω) = 0", lw = 2, color = :blue)
plot!(verLine, LinRange(minimum(y), maximum(y), length(y)), label = "x = $c", lw = 2, color = :blue)
scatter!([c], [0], label = "Root", color = :blue)