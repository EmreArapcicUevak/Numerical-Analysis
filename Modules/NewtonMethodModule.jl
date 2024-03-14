module NewtonMethodModule
        using LinearAlgebra
        export NewtonMethod, QuasiNewtonMethod, MultiDimentionalNewtonMethod, AproximateJacobian
    
        function NewtonMethod(f :: Function, f_prime :: Function, x₀ :: Real; δ = 1e-15, ϵ = 1e-15, maxIterations = 1000)
            for i ∈ 1:maxIterations
                local f_prime_x₀, f_x₀ = f_prime(x₀), f(x₀)

                #=
                if abs(f_prime_x₀) < ϵ
                    error("Derivative is too small, indicating a possible horizontal tangent or numerical instability at x = $x₀")
                end
                =#

                local h = f_x₀ / f_prime_x₀
                if abs(h / x₀) ≤ ϵ || abs(f_x₀) ≤ δ
                    return (c = x₀,iterations = i)
                end

                x₀ -= h 
            end

            return nothing
        end

        function QuasiNewtonMethod(f :: Function, x₀ :: Real, x₁ :: Real; δ = 1e-15, ϵ = 1e-15, maxIterations = 1000)
            local f₀, f₁ = f(x₀), f(x₁)
            
            for i ∈ 1:maxIterations
                local xₖ = x₁ - (x₁ - x₀) * f₁ / (f₁ - f₀)
                x₀, x₁ = x₁, xₖ
                f₀, f₁ = f(x₀), f(x₁)

                println("k: $i, xₖ: $xₖ, f(xₖ): $(f(xₖ))")
                if abs(x₁ - x₀) ≤ ϵ || abs(f(x₁)) ≤ δ
                    return (c = x₁,iterations = i)
                end
            end

            return nothing
            
        end

        function MultiDimentionalNewtonMethod(F :: Function, J :: Function, x₀ :: Vector{T}; δ = 1e-15, ϵ = 1e-15, maxIterations = 1000) where T <: Real
            for i ∈ 1:maxIterations
                local J_x₀, F_x₀ = J(x₀), F(x₀)
                local x₁ = x₀ - J_x₀ \ F_x₀

                if norm(x₁ - x₀, 2) ≤ ϵ || norm(F(x₁), 2) ≤ δ
                    return (c = x₁,iterations = i)
                end

                x₀ = x₁
            end

            return nothing
            
        end

        function AproximateJacobian(F :: Function, x₀ :: Vector{T}; t = 1e-6 :: Float64) where T <: Real
            local N = length(x₀)
            local J = zeros(T, N, N)

            for i ∈ 1:N
                local x = copy(x₀)
                x[i] += t
                J[:,i] = (F(x) - F(x₀)) / t
            end

            return J
        end
end