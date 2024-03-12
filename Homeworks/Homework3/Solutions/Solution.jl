push!(LOAD_PATH, pwd()*"/Modules/")
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

function F(X :: Vector)
    N = length(X)
    Fₓ = Vector{Float64}()
    h = 1/(N+1)
    push!(Fₓ, (X[2] - 2X[1])/h^2 + exp(X[1]))
    for i ∈ 2:N-1
        push!(Fₓ, (X[i-1] + X[i+1] - 2*X[i])/h^2 + exp(X[i]))
    end
    push!(Fₓ, (X[N-1] - 2X[N])/h^2 + exp(X[N]))

    return Fₓ
end

function J(X :: Vector)
    N = length(X)
    Jₓ, h = zeros(Float64,N, N), 1/(N+1)
    c = 1/h^2

    Jₓ[1,1] = -2/h^2 + exp(X[1])
    Jₓ[1,2] = c

    Jₓ[N,N] = -2/h^2 + exp(X[N])
    Jₓ[N,N-1] = c

    for i ∈ 2:N-1
        Jₓ[i,i] = -2/h^2 + exp(X[i])
        Jₓ[i,i-1] = c
        Jₓ[i,i+1] = c
    end

    return Jₓ
end


for N ∈ 140:10:200
    c, itt = MultiDimentionalNewtonMethod(F, J, zeros(N); maxIterations = 9*10^6)
    scatter(LinRange(0,1,N+2), [0; c; 0], label = "N = $N", ylabel = "y(t)", xlabel = "t")
    savefig("Plot$(N).png")
    println("N = $N, iterations = $itt")
end