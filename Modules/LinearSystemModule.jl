module LinearSystemModule

    using MatrixAnalysisModule, ElementaryMatrixOperationsModule
    export SolveUpperDiagonal, SolveLowerDiagonal

    validType = Union{Rational, AbstractFloat}

    function SolveUpperDiagonal(A :: Matrix{T}, b :: Vector{T}) where T <: validType
        @assert IsUpperTriangular(A) "A must be upper triangular"

        n = length(b)
        x = zeros(T, n)
        x[n] = b[n] / A[n,n]

        for i ∈ n-1:-1:1
            x[i] = (b[i] -  sum([A[i,j] * x[j] for j ∈ i+1:n]) ) / A[i,i]
        end

        return x
    end

    function SolveLowerDiagonal(A :: Matrix{T}, b :: Vector{T}) where T <: validType
        @assert IsLowerTriangular(A) "A must be lower triangular"

        n = length(b)
        x = zeros(T, n)
        x[1] = b[1] / A[1,1]

        for i ∈ 2:n
            x[i] = (b[i] -  sum([A[i,j] * x[j] for j ∈ 1:i-1]) ) / A[i,i]
        end

        return x
    end

    export GetIdentityMatrix
    function GetIdentityMatrix(T :: Type,n)
        I = zeros(T, n, n)
        for i ∈ 1:n
            I[i,i] = 1
        end
        
        return I
    end

    export PLUFactorization!
    function PLUFactorization!(A :: Matrix{T}) :: NamedTuple{(:P, :L), Tuple{Vector{Int64}, Matrix{T}}} where T <: validType
        rows, columns = size(A)
        L = GetIdentityMatrix(T, rows)
        P = collect(1:columns)

        for i ∈ 1:columns
            swapIndex = argmax(abs.(A[i:end,i])) + i - 1
            if i != swapIndex
                SwapRows!(A, i, swapIndex)
                L[i,1:i-1], L[swapIndex,1:i-1] = L[swapIndex,1:i-1], L[i,1:i-1]
                P[i], P[swapIndex] = P[swapIndex], P[i]
            end

            @assert A[i,i] != 0 "A is singular"

            for j ∈ i+1:rows
                c = A[j,i] / A[i,i]
                AddRowMultiple!(A, j, i, -c)
                L[j,i] = c
            end
        end
        
        return (P = P, L = L)
    end

    export PLUFactorization
    function PLUFactorization(A :: Matrix{T}) :: NamedTuple{(:P, :U, :L), Tuple{Vector{Int64}, Matrix{T}, Matrix{T}}}  where T <: validType
        local U = copy(A)
        local P, L = PLUFactorization!(U)
        return (P = P, U = U, L = L)
    end

    export SolvePLU
    function SolvePLU(P :: Vector{Int64}, L :: Matrix{T}, U :: Matrix{T}, b :: Vector{T}) where T <: validType
        y = SolveLowerDiagonal(L, b[P])
        x = SolveUpperDiagonal(U, y)
        
        return x
    end

    function SolvePLU(Result :: NamedTuple{(:P, :U, :L), Tuple{Vector{Int64}, Matrix{T}, Matrix{T}}}, b :: Vector{T}) where T <: validType 
        return SolvePLU(Result.P, Result.L, Result.U, b)
    end

    export norm
    function norm(x :: Vector{T}, p :: Int) where T <: validType
        @assert p > 0 "p must be greater than 0"
        return sum(abs.(x).^p)^(1/p)
    end
end