push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

#### Other modules that you might need go here ####
using NewtonMethodModule
using LinearSystemModule
###################################################



#### Your code goes here ####
N = 200
const α, k, m , l₀, xd, a = 1, 1, 1, √2,  Float64[0, 2],  Float64[0, 1]
const x₀, v₀, λₙ, μₙ = Float64[1, 1],  Float64[0, 0],  Float64[0, 0],  Float64[0, 0]
Δt = 10 / N

function Lᵢ(xᵢ :: Vector{T}, uᵢ :: T) :: T where T <: Real
    norm(xᵢ - [uᵢ,0], 2)
end

function F(X :: Vector{T}) :: Vector{T} where T <: Real
    local N = (length(X) - 1) ÷ 9

    local result = Float64[]
    
    local x = reshape(X[1:2N], 2, N)
    local v = reshape(X[2N + 1:4N], 2, N)
    local λ = reshape(X[4N + 1:6N], 2, N)
    local μ = reshape(X[6N + 1:8N], 2, N)
    local u = X[8N + 1:9N + 1]

    for i ∈ 1:N-2
        local r = (λ[:, i+2] - λ[:, i+1]) / Δt + μ[:, i+2] # Eq 1
        push!(result, r[1])
        push!(result, r[2])

        r = begin # Eq 2
            local d = (x[:, i] - [u[i+1], 0])
            local Lₙ = Lᵢ(x[:, i], u[i+1])
            local c₁ = (μ[:, i+2] - μ[:, i+1]) / Δt
            local c₂ = x[:, i] - xd
            local c₃ = k * l₀ * d * d' / (m * Lₙ^3)
            local c₄ = k/m * (Lₙ - l₀)/Lₙ * [1 0; 0 1]
            c₁ - c₂ - (c₃ + c₄) * λ[:, i+2]
        end
        push!(result, r[1])
        push!(result, r[2])

        r = begin # Eq 3
            local d = (x[:, i] - [u[i+1], 0])
            local Lₙ = Lᵢ(x[:, i], u[i+1])
            local c₁ = (v[:, i+1] - v[:, i]) / Δt
            local c₂ = k * (Lₙ - l₀) * d / (m * Lₙ)
            local c₃ = a / m
            c₁ + c₂ - c₃
        end
        push!(result, r[1])
        push!(result, r[2])

        local r = begin # Eq 4
            local c₁ = (x[:, i+1] - x[:, i]) / Δt
            c₁ - v[:, i]
        end

        push!(result, r[1])
        push!(result, r[2])

        r = begin # Eq 5
            local d = (x[:, i] - [u[i+1], 0])
            local Lₙ = Lᵢ(x[:, i], u[i+1])

            local c₁ = α*u[i+1]
            local c₂ = (-k * l₀ * [1 0] * d / (m * Lₙ^3))[1] * d
            local c₃ = -k / m * (Lₙ - l₀) / Lₙ * [1, 0]
            c₁ + λ[:, i+1]' * (c₂ + c₃)
        end

        push!(result, r)
    end

    ############################## First equation ##############################
    local rSpecial = (λ[:, 2] - λ[:, 1]) / Δt + μ[:, 2]
    push!(result, rSpecial[1])
    push!(result, rSpecial[2])

    rSpecial = (λₙ - λ[:, N]) / Δt + μₙ
    push!(result, rSpecial[1])
    push!(result, rSpecial[2])
    #############################################################################

    ############################## Second equation ##############################
    rSpecial = begin
        local d = (x₀ - [u[1], 0])
        local Lₙ = Lᵢ(x₀, u[1])
        local c₁ = (μ[:, 2] - μ[:, 1]) / Δt
        local c₂ = x₀ - xd
        local c₃ = k * l₀ * d * d' / (m * Lₙ^3)
        local c₄ = k/m * (Lₙ - l₀)/Lₙ * [1 0; 0 1]
        c₁ - c₂ - (c₃ + c₄) * λ[:, 2]
    end
    push!(result, rSpecial[1])
    push!(result, rSpecial[2])

    rSpecial = begin
        local d = (x[:, N-1] - [u[N], 0])
        local Lₙ = Lᵢ(x[:, N-1], u[N])
        local c₁ = (μₙ - μ[:, N]) / Δt
        local c₂ = x[:, N-1] - xd
        local c₃ = k * l₀ * d * d' / (m * Lₙ^3)
        local c₄ = k/m * (Lₙ - l₀)/Lₙ * [1 0; 0 1]
        c₁ - c₂ - (c₃ + c₄) * λₙ
    end
    push!(result, rSpecial[1])
    push!(result, rSpecial[2])
    #############################################################################

    ############################## Third equation ##############################
    rSpecial = begin
        local d = (x₀ - [u[1], 0])
        local Lₙ = Lᵢ(x₀, u[1])
        local c₁ = (v[:, 1] - v₀) / Δt
        local c₂ = k * (Lₙ - l₀) * d / (m * Lₙ)
        local c₃ = a / m
        c₁ + c₂ - c₃
    end
    push!(result, rSpecial[1])
    push!(result, rSpecial[2])

    rSpecial = begin
        local d = (x[:, N-1] - [u[N], 0])
        local Lₙ = Lᵢ(x[:, N-1], u[N])
        local c₁ = (v[:, N] - v[:, N-1]) / Δt
        local c₂ = k * (Lₙ - l₀) * d / (m * Lₙ)
        local c₃ = a / m
        c₁ + c₂ - c₃
    end

    push!(result, rSpecial[1])
    push!(result, rSpecial[2])
    #############################################################################

    ############################## Fourth equation ##############################
    rSpecial = begin
        local c₁ = (x[:, 1] - x₀) / Δt
        c₁ - v₀
    end
    push!(result, rSpecial[1])
    push!(result, rSpecial[2])

    rSpecial = begin
        local c₁ = (x[:, N] - x[:, N-1]) / Δt
        c₁ - v[:, N-1]
    end

    push!(result, rSpecial[1])
    push!(result, rSpecial[2])

    ############################ Fifth Equation  #############################
    rSpecial = begin
        local d = (x₀ - [u[1], 0])
        local Lₙ = Lᵢ(x₀, u[1])
        local c₁ = α*u[1]
        local c₂ = (-k * l₀ * [1 0] * d / (m * Lₙ^3))[1] * d
        local c₃ = -k / m * (Lₙ - l₀) / Lₙ * [1, 0]
        c₁ + λ[:, 1]' * (c₂ + c₃)
    end

    push!(result, rSpecial)

    rSpecial = begin
        local d = (x[:, N-1] - [u[N], 0])
        local Lₙ = Lᵢ(x[:, N-1], u[N])

        local c₁ = α*u[N]
        local c₂ = (-k * l₀ * [1 0] * d / (m * Lₙ^3))[1] * d
        local c₃ = -k / m * (Lₙ - l₀) / Lₙ * [1, 0]
        c₁ + λ[:, N]' * (c₂ + c₃)
    end

    push!(result, rSpecial)

    rSpecial = begin
        local d = (x[:, N] - [u[N + 1], 0])
        local Lₙ = Lᵢ(x[:, N], u[N + 1])
        local c₁ = α*u[N + 1]
        local c₂ = (-k * l₀ * [1 0] * d / (m * Lₙ^3))[1] * d
        local c₃ = -k / m * (Lₙ - l₀) / Lₙ * [1, 0]
        c₁ + λₙ' * (c₂ + c₃)
    end

    push!(result, rSpecial)
    ##########################################################################

    @assert length(result) == 9N + 1 "Size missmatch"
    return result
end

guess = rand(Float64,9N + 1);

J(X :: Vector) = AproximateJacobian(F, X)
result = MultiDimentionalNewtonMethod(F, J, guess; δ = 5e-10)
u = result.c[8N + 1:end];
x = LinRange(0,10, length(u));
plot(x, u, label = "u(t)", title = "Optimal Control α = $(α)", xlabel = "t", ylabel = "u(t)", line = 2)
savefig("$(N)Points$(α)Alpha.png")