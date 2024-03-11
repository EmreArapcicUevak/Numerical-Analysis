push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

function F(X :: Vector)
    N = length(X)
    Fₓ = Vector{Float64}()
    h = 1/(N+2)
    push!(Fₓ, (X[2] - 2X[1])/h^2 + exp(X[1]))
    for i ∈ 2:N-1
        push!(Fₓ, (X[i-1] + X[i+1] - 2*X[i])/h^2 + exp(X[i]))
    end
    push!(Fₓ, (X[N-2] - 2X[N])/h^2 + exp(X[N]))

    return Fₓ
end

function J(X :: Vector)
    N = length(X)
    Jₓ, h = zeros(Float64,N, N), 1/(N+2)
    for i ∈ 1:N
        Jₓ[i,i] = -2/h^2 + exp(X[i])
        if i > 1
            Jₓ[i, i-1] = 1/h^2
        end
        if i < N
            Jₓ[i, i+1] = 1/h^2
        end
    end

    return Jₓ
end

N = 100
c, itt = MultiDimentionalNewtonMethod(F, J, zeros(N); maxIterations = 9*10^6)

scatter(LinRange(0,1,N+2), [0; c; 0], label = "y(t)", ylabel = "y(t)", xlabel = "t")

# Begin LaTeX table code
latex_code = """
\\begin{table}[ht]
\\centering
\\begin{tabular}{|c|c|}
\\hline
Index & Value \\\\
\\hline
"""

# Populate table with y_values
for (i, y) in enumerate([0;c;0])
    latex_code *= "$i & $y \\\\\n\\hline\n"
end

# End LaTeX table code
latex_code *= """
\\end{tabular}
\\caption{Table of Points}
\\end{table}
"""

# Output LaTeX code (you can redirect this to a .tex file or print it)
println(latex_code)