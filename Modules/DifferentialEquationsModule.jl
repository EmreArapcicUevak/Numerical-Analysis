module DifferentialEquationsModule
    using NewtonMethodModule

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

    function ExplicitEulerMethod(f :: Function, z₀ :: Vector, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Vector{Real}}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n
        z = Vector{Real}[z₀]

        for i ∈ 1:n
            push!(z, z[end] + h * f(bounds[begin] + (i - 1)*h, z[end]))
        end

        return z
    end

    export ImplicitEulerMethod
    function ImplicitEulerMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n
        y = Real[y₀]
        
        for i ∈ 1:n
            g(x :: Real) = x - y[end] - h * f(bounds[begin] + i*h, x)
            y₀ = QuasiNewtonMethod(g, y[end], 1e-6)
                
            @assert y₀ !== nothing "Newton's method failed to converge"
            push!(y, y₀)
        end

        return y
    end

    export ModifiedEulerMethod
    function ModifiedEulerMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n

        y = Real[y₀]
        for i ∈ 1:n
            k₁ = f(bounds[begin] + (i - 1)*h, y[end])
            y₁ = y[end] + h*k₁

            k₂ = f(bounds[begin] + i*h, y₁)
            push!(y, y[end] + h/2 * (k₁ + k₂))
        end

        return y
    end

    function ModifiedEulerMethod(f :: Function, z₀ :: Vector, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Vector{Real}}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n
        z = Vector{Real}[z₀]

        for i ∈ 1:n
            k₁ = f(bounds[begin] + (i - 1)*h, z[end])
            z₁ = z[end] + h*k₁

            k₂ = f(bounds[begin] + i*h, z₁)
            push!(z, z[end] + h/2 * (k₁ + k₂))
        end

        return z
        
    end

    export MidpointMethod
    function MidpointMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n

        y = Real[y₀]
        for i ∈ 1:n
            t₀ = bounds[begin] + (i - 1)*h
            tₘ = t₀ + h/2

            yₘ = y[end] + h/2 * f(t₀, y[end])
            y₁ = y[end] + h * f(tₘ, yₘ)

            push!(y, y₁)
        end

        return y
    end

    function MidpointMethod(f :: Function, z₀ :: Vector, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Vector{Real}}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n
        z = Vector{Real}[z₀]

        for i ∈ 1:n
            t₀ = bounds[begin] + (i - 1)*h
            tₘ = t₀ + h/2

            zₘ = z[end] + h/2 * f(t₀, z[end])
            z₁ = z[end] + h * f(tₘ, zₘ)

            push!(z, z₁)
        end

        return z
        
    end

    #=
    export RungeKuttaMethod
    function RungeKuttaMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n

        y = Real[y₀]
        for i ∈ 1:n
            t₀ = bounds[begin] + (i - 1)*h
            k₁ = f(t₀, y[end])
            k₂ = f(t₀ + h/2, y[end] + h/2 * k₁)
            k₃ = f(t₀ + h/2, y[end] + h/2 * k₂)
            k₄ = f(t₀ + h, y[end] + h * k₃)

            push!(y, y[end] + h/6 * (k₁ + 2*k₂ + 2*k₃ + k₄))
        end

        return y
    end

    export GeneralizedRungeKuttaMethod
    function GeneralizedRungeKuttaMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}, c1 :: Real, c2 :: Real) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        @assert c1 + c2 == 1 "c₁ + c₂ has to be equal to 1"
        @assert c1 != 0 && c2 != 0 "c₁ and c₂ cannot be equal to 0"
        
        # c1, c2 ∈ ℝ, α, β = 1/(2*c2), 1/(2*c1)
        α, β = 1/(2*c2), 1/(2*c1)
        h = (bounds[2] - bounds[1]) / n

        y = Real[y₀]
        for i ∈ 1:n

        end

        return y
    end

    export HeunMethod
    function HeunMethod(f :: Function, y₀ :: Real, n :: Int, bounds :: Tuple{Real, Real}) :: Vector{Real}
        @assert n > 0 "N has to be greater than 0"
        @assert bounds[1] < bounds[2] "a has to be less than b"
        h = (bounds[2] - bounds[1]) / n

        y = Real[y₀]
        for i ∈ 1:n
            t₀ = bounds[begin] + (i - 1)*h
            k₁ = f(t₀, y[end])
            k₂ = f(t₀ + h, y[end] + h*k₁)

            push!(y, y[end] + h/2 * (k₁ + k₂))
        end

        return y
    end
    =#
end