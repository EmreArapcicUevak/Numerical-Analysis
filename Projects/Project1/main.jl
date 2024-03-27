push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, CairoMakie# Always leave this on!

#### Other modules that you might need go here ####
using NewtonMethodModule
###################################################

#### Your code goes here ####
const startingPoint, endingPoint = 1,2
const yStart, yEnd = 0, log(2)

const numericalTypes = Union{AbstractFloat, Rational}

function F(Y :: Vector{T}) :: Vector{T} where T <: numericalTypes
    local X, N = T[], length(Y)
    @assert N ≥ 2 "Y must have at least 2 points"

    local h = (endingPoint - startingPoint) / (N + 1)
    local h² = h^2

    push!(X, begin
        local d²y = (Y[2] - 2Y[1] + yStart) / h²
        local dy = (Y[2] - yStart) / (2h)
        local tᵢ = startingPoint + h
        d²y - dy^2 + Y[1] - log(tᵢ)
    end)

    for i ∈ 2:N-1
        push!(X, begin
            local d²y = (Y[i+1] - 2Y[i] + Y[i-1]) / h²
            local dy = (Y[i+1] - Y[i-1]) / (2h)
            local tᵢ = startingPoint + i*h
            d²y - dy^2 + Y[i] - log(tᵢ)
        end)
    end

    push!(X, begin
        local d²y = (yEnd - 2Y[N] + Y[N-1]) / h²
        local dy = (yEnd - Y[N-1]) / (2h)
        local tᵢ = endingPoint - h
        d²y - dy^2 + Y[N] - log(tᵢ)
    end)

    return X
end


function J(X :: Vector{T}) :: Matrix{T} where T <: numericalTypes
    return AproximateJacobian(F, X)
end

figure = Figure(backgroundcolor = "#c8d6e5", size=(900, 1500));
for (index, N) ∈ enumerate(collect(10:10:400))
    t = LinRange(startingPoint, endingPoint, N + 2)
    solution = MultiDimentionalNewtonMethod(F, J, zeros(Float64, N))

    if solution !== nothing
        yPoints = [yStart; solution.c; yEnd]
        row = div(index - 1, 4) + 1
        column = mod(index - 1, 4) + 1
        ax = Axis(figure[row, column], title = "N = $N", xlabel = "t", ylabel = "y(t)")
        lines!(ax, t, yPoints; color = "#54a0ff")
    else
        println("No solution found for N = $N")
        break
    end
end

figure
save("Projects/Project1/NewtonMethod.png", figure, px_per_unit = 2)