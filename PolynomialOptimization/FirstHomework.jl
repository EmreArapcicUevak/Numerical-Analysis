function f(n :: Integer, coefficients :: Union{AbstractArray, Tuple}, x :: Real)
    @assert length(coefficients) ≥ n
    local value = 0

    for (index, coef) ∈ enumerate(coefficients)
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
for n ∈ map(x -> 10^x, 2:9)
    println("Function comparisson for a $n degree polynomial")
    
    for x₀ ∈ LinRange(0,1,3)
        local result, time, garbage... = @timed fOptimized(n, randomCoefficients, x₀) 
        println("f₁($(x₀)) = $(result): Time taken $(time) seconds")

        local result, time, garbage... = @timed f(n, randomCoefficients, x₀) 
        println("f₂($(x₀)) = $(result): Time taken $(time) seconds")
        println("-"^20)
    end

    println("*"^50)
end 