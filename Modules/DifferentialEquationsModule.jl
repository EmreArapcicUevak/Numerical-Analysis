module DifferentialEquationsModule
    export  ExplicitEulerMethod
    function ExplicitEulerMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n
        y = Real[y₀]
        
        for i ∈ 1:n
            push!(y, y[end] + h * f(bounds[begin] + (i - 1)*h, y[end]))
        end

        return y
    end
end