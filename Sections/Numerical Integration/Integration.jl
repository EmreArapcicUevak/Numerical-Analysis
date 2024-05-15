push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

using NumericalIntegrationModule

function PlotAbsoluteError(funcToIntegrate :: Function, a :: Real, b :: Real, exactValue :: Real, n :: StepRange{Int64, Int64}, rules::Array{Function, 1}, labels::Array{String, 1}; activationFunction :: Function = log10)
    @assert length(rules) == length(labels) "The number of rules and labels must be the same"

    local x = collect(n)
    local firstPlot = true

    for (rule, label) ∈ zip(rules, labels)
        local Error = Real[]
        for i ∈ n
            local absError = abs(rule(funcToIntegrate, a, b, n = i) - exactValue)
            push!(Error, absError)
        end
        
        if firstPlot 
            firstPlot = false
            scatter(x, activationFunction.(Error), label = label, xlabel = "n", ylabel = "Absolute Error", title = "Absolute Error of rules", lw = 2)
        else
            scatter!(x, activationFunction.(Error), label = label, lw = 2)
        end
    end

    display(plot!())
end

f(x) = sin(x)
PlotAbsoluteError(f, 0, pi, 2, 10:10:1000, [LeftRectangleRule, MidpointRule, TrapezoidalRule, SimpsonRule], ["Left Rectangle Rule", "Midpoint Rule", "Trapezoidal Rule", "Simpson Rule"]; activationFunction = x -> x)


f(x) = exp(x)  
PlotAbsoluteError(f, 0, 1, exp(1) - 1, 3:3:3*100, [LeftRectangleRule, MidpointRule, TrapezoidalRule, SimpsonRule3_8], ["Left Rectangle Rule", "Midpoint Rule", "Trapezoidal Rule", "Simpson Rule"])


f(x) = x^3
PlotAbsoluteError(f, 0, 1, 1/3, 9:3:3 * 100, [LeftRectangleRule, MidpointRule, TrapezoidalRule, SimpsonRule3_8], ["Left Rectangle Rule", "Midpoint Rule", "Trapezoidal Rule", "Simpson Rule 3/8"])

f(x) = sin(x)
PlotAbsoluteError(f, 0, pi, 2, 9:3:3 * 100, [LeftRectangleRule, MidpointRule, TrapezoidalRule, SimpsonRule3_8], ["Left Rectangle Rule", "Midpoint Rule", "Trapezoidal Rule", "Simpson Rule"])

savefig("AbsoluteError.png")