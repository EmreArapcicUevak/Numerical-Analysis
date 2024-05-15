push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
Plots.plotlyjs() # Set PlotlyJS as the backend for Plots.jl

#### Other modules that you might need go here ####
using NewtonMethodModule
###################################################

#### Your code goes here ####

function SolveLinearODE(a₁ :: Function, a₀ :: Function, f :: Function;boundaries :: Tuple{Float64, Float64} ,range :: Tuple{Float64, Float64}, N :: Int)
    @assert N > 0 "The number of points must be greater than 0"
    local α, β = boundaries
    local a,b = range
    local h = (b - a) / (N + 1)

    function VectorField(Y :: Vector)
        Y = [α;Y;β]    
        
        local Yₙ = zeros(N)
        for i ∈ 1:N
            Yₙ[i] = begin
                y⁽²⁾ = (Y[i + 2] - 2Y[i + 1] + Y[i]) / h^2
                y⁽¹⁾ = (Y[i + 2] - Y[i]) / (2h)
                xᵢ = a + i * h
                y⁽²⁾ + a₁(xᵢ) * y⁽¹⁾ + a₀(xᵢ) * Y[i + 1] - f(xᵢ)
            end
        end

        return Yₙ
    end

    return MultiDimentionalNewtonMethod(VectorField, x -> AproximateJacobian(VectorField, x), zeros(N), δ = 1e-15, ϵ = 1e-15, maxIterations = 1000)
end

a₁(x) = sin(x^2) + exp(x)
a₀(x) = cos(x^5) + log(x)
f(x) = 2x^2 - 2x + 1

boundaries = (.5, 1.) ## VALUE OF THE FUNCTION AT THE BOUNDARIES
range = (0., 2.) ## RANGE OF THE FUNCTION
N = 100

result = SolveLinearODE(a₁, a₀, f, boundaries = boundaries, range = range, N = N)

plot(LinRange(range[1], range[2], N+2), [boundaries[1] ;result.c ; boundaries[2]], label = "Numerical Solution", xlabel = "x", ylabel = "y", title = "Numerical Solution of the ODE", line = :dash, lw = 2, color = :red, legend = :topleft)
scatter!(LinRange(range[1], range[2], N+2), [boundaries[1] ;result.c ; boundaries[2]], label = "Numerical Solution", color = :blue, markersize = 2)