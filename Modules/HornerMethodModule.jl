module HornerMethodModule
    export hornerMethodPolynomial

    function hornerMethodPolynomial(p :: Vector{T}, x :: Real) where T <: Real
        local n = length(p)
        local result = p[n]
        for i ∈ n-1:-1:1
            result = result*x + p[i]
        end
        return result
    end

end