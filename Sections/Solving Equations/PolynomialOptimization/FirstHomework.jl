using Plots
function f(n :: Integer, coefficients :: Union{AbstractArray, Tuple}, x :: Real)
    @assert length(coefficients) ≥ n
    local value = 0

    for (index, coef) ∈ enumerate(coefficients[begin:n])
        value += coef * x ^ (index - 1)
    end

    return value 
end


function fOptimized(n :: Integer, coefficients :: Union{AbstractArray, Tuple}, x :: Real)
    @assert length(coefficients) ≥ n

    local value = coefficients[n]

    for i ∈ n-1:-1:1
        value = coefficients[i] + value*x
    end

    return value
end

randomCoefficients = rand(1_000_000_000)

allDegrees = map(x -> 10^x, 1:9) # You can change this for a smaller range if you want to test the code faster or add more points for better scatters.
time = Dict(
    "Optimized" => [],
    "Unoptimized" => []
)

for (index, n) ∈ enumerate(allDegrees)
    println("Function comparisson for a $n degree polynomial")
    
    local totalTime = [0.0, 0.0]

    for x₀ ∈ 0:0.5:1
        local result, time, garbage... = @timed fOptimized(n, randomCoefficients, x₀) 
        totalTime[1] = time > totalTime[1] ? time : totalTime[1]

        println("f₁($(x₀)) = $(result): Time taken $(time) seconds")

        local result, time, garbage... = @timed f(n, randomCoefficients, x₀) 
        totalTime[2] = time > totalTime[2] ? time : totalTime[2]

        println("f₂($(x₀)) = $(result): Time taken $(time) seconds")
        println("-"^20)
    end

    append!(time["Optimized"], totalTime[1])
    append!(time["Unoptimized"], totalTime[2])
    println("*"^20)
end 

scatter(allDegrees, time["Optimized"], label="Optimized", xlabel="Degree of the polynomial", ylabel="Time taken (s)", title="Optimized vs Unoptimized function")
scatter!(allDegrees, time["Unoptimized"], label="Unoptimized")