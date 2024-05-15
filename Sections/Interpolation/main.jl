push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots # Always leave this on!
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

#### Other modules that you might need go here ####

###################################################
using LinearSystemModule
using InterpolationModule
#### Your code goes here ####

A = [1 0 0 0; 1 π/2 (π/2)^2 (π/2)^3; 1 π (π)^2 (π)^3; 1 3π/2 (3π/2)^2 (3π/2)^3]
b = Float64[0, 1, 0, -1]

result = PLUFactorization(A)
coef = SolvePLU(result, b)

x = LinRange(0, 6, 9000)
y = sin.(x)
poly = coef[1] .+ coef[2] .* x .+ coef[3] .* x.^2 .+ coef[4] .* x.^3

plot(x, y, label = "sin(x)")
scatter!([0, pi/2, pi, 3pi/2],[0,1,0,-1], label = "Value Points")
plot!(x, poly, label = "Interpolated Polynomial")

x = Float64[2, 2.75, 4]
y = 1 ./ x
polyCoef = RecursiveInterpolation(x, y)
poly(X) = polyCoef[1] + polyCoef[2] * (X - x[1]) + polyCoef[3] * (X - x[1]) * (X - x[2])
poly(3)
poly(2)
A = [1 2 4; 1 2.75 2.75^2; 1 4 4^2]
result = PLUFactorization(A)
coef = SolvePLU(result, y)
coef

polyGuess = LagrangeInterpolation(x, y)
polyGuess(3)

coef = LinearInterpolation(x, y)

x = Float64[0,1,2,3]
y = Float64[3,5,7,9]
coef = NewtonDivider(x, y)

poly(X) = coef[1] + coef[2] * (X - x[1]) + coef[3] * (X - x[1]) * (X - x[2]) + coef[4] * (X - x[1]) * (X - x[2]) * (X - x[3])
xExp = LinRange(0, 3, 9000)
yExp = poly.(xExp)

scatter(x, y, label = "Value Points")
plot!(xExp, yExp, label = "Interpolated Polynomial")


## Chebychev Interpolation
f(x :: Real) = 1 / (1 + 25 * x^2)
N = 15
a,b = -pi/4,pi/4
x = LinRange(a,b , 10000)
plot(x, f.(x), label = "Value Points")

xᵢ = ChebyshevNodes(a,b, N)
coefCheb = NewtonDividedDifference(xᵢ, f.(xᵢ))
polyAprox = NewtonInterpolation(xᵢ, coefCheb)

plot!(x, polyAprox.(x), label = "Chebychev Polynomial $(N)")

xⱼ = collect(LinRange(a,b, N))

coefLinDistance = NewtonDividedDifference(xⱼ, f.(xⱼ))
polyLinDistance = NewtonInterpolation(xⱼ, coefLinDistance)
plot!(x, polyLinDistance.(x), label = "Linearly Spaced Polynomial $(N)")

## Hermitian Interpolation
f(x :: Real) = exp(x) + sin(10x)
fPrime(x :: Real) = exp(x) + 10cos(10x)

x = Float64[0, 0.4, 1, 2, 2.6, 3]
y = f.(x)
yPrime = fPrime.(x)

hermitCoef = HermitianInterpolationCoef(x, y, yPrime)

xᵢ = LinRange(0, 3, 100)
yᵢ = f.(xᵢ)
InterpoaltionGuess = NewtonInterpolation([X for X ∈ x for _ ∈ 1:2], hermitCoef)

plot(xᵢ, yᵢ, label = "Value Points")
plot!(xᵢ, InterpoaltionGuess.(xᵢ), label = "Hermitian Interpolation")
InterpoaltionGuess.(x)
y
