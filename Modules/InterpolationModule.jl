module InterpolationModule
    using LinearSystemModule

    export RecursiveInterpolation
    function RecursiveInterpolation(X :: Vector{T1}, Y :: Vector{T2}) :: Vector{Float64} where {T1 <: Real, T2 <: Real}
        local N = length(X)
        @assert N == length(Y) "X and Y must have the same length"
        
        local CMatrix = zeros(Float64, N, N)
        CMatrix[:,1] .= 1

        for row ∈ 2:N
            for column ∈ 2:row
                CMatrix[row, column] = CMatrix[row, column - 1] * (X[row] - X[column - 1])
            end
        end

        println(CMatrix)
        return SolveLowerDiagonal(CMatrix, Y)
    end

    export LagrangeInterpolation
    function LagrangeInterpolation(X :: Vector{T}, Y :: Vector{T}) :: Function where T <: Real
        local N = length(X)
        @assert N == length(Y) "X and Y must have the same length"

        return x :: Real -> sum([Y[i] * prod([(x - X[j]) / (X[i] - X[j]) for j ∈ 1:N if j != i]) for i ∈ 1:N])
    end

    export LinearInterpolation
    function LinearInterpolation(X :: Vector{T}, Y :: Vector{T}) where T <: Real
        local N = length(X)
        @assert N == length(Y) "X and Y must have the same length"

        A = zeros(Float64, N, N)
        A[:,1] .= 1

        for i ∈ 2:N
            for j ∈ 2:N
                A[i, j] = A[i, j-1] * X[i]
            end
        end

        return SolvePLU(PLUFactorization(A), Y)
    end

    export NewtonDividedDifference!
    function NewtonDividedDifference!(X :: Vector{T}, Y :: Vector{T}) where T <: AbstractFloat
        local N = length(X)
        @assert N == length(Y) "X and Y must have the same length"

        for i ∈ 1:N
            for j ∈ N:-1:i+1
                Y[j] = (Y[j] - Y[j-1]) / (X[j] - X[j-i])
            end
        end
    end

    export NewtonDividedDifference
    function NewtonDividedDifference(X :: Vector{T}, Y :: Vector{T}) :: Vector{T} where T <: AbstractFloat
        local yCopy = copy(Y)
        NewtonDividedDifference!(X, yCopy)
        
        return yCopy
    end

    export NewtonInterpolation
    function NewtonInterpolation(X :: Vector{T}, Coef :: Vector{T}) :: Function where T <: Real
        local N = length(X)
        @assert N == length(Coef) "X and Coef must have the same length"


        return x :: Real -> sum([Coef[i] * prod([x - X[j] for j ∈ 1:i-1]) for i ∈ 1:N])
    end

    export ChebyshevNodes
    function ChebyshevNodes(a :: Real, b :: Real, N :: Int) :: Vector{Float64}
        return [(b + a) * 0.5 + (b - a) * 0.5 * cos((2i - 1) * π / (2N)) for i ∈ 1:N]
    end
end