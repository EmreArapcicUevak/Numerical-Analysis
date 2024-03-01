module NewtonMethodModule
        export NewtonMethod, QuasiNewtonMethod
    
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

                println("k: $i, xₖ: $xₖ, f(xₖ): $(f(xₖ))")
                if abs(x₁ - x₀) ≤ ϵ || abs(f(x₁)) ≤ δ
                    return (c = x₁,iterations = i)
                end
            end

            return nothing
            
        end
end