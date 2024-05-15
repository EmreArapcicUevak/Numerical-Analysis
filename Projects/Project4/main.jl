push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
Plots.plotlyjs() # Set PlotlyJS as the backend for Plots.jl

#### Other modules that you might need go here ####
using NumericalIntegrationModule
###################################################

#### Your code goes here ####

function DoubleIntegral(f :: Function, firstBoundaries :: Tuple{Real,Real}, secondBoundaries :: Tuple{Function, Function};IntegrationRule :: Function = TrapezoidalRule , n :: Int = 1000)
    local (α, β) = secondBoundaries
    local (a, b) = firstBoundaries
    function g(x)
        LowerBound, UpperBound = α(x), β(x)

        local modifiedFunc(y) = f(x, y)
        return IntegrationRule(modifiedFunc, LowerBound, UpperBound; n = n)
    end

    return IntegrationRule(g, a, b; n = n)
end

f(r, θ) = r
xBoundaries = (0, 1)
yBoundaries = (x -> 0, x -> 2*pi)
DoubleIntegral(f, xBoundaries, yBoundaries; n = 300)
DoubleIntegral(f, xBoundaries, yBoundaries;IntegrationRule = SimpsonRule, n = 300)
DoubleIntegral(f, xBoundaries, yBoundaries;IntegrationRule = SimpsonRule3_8, n = 300)
π