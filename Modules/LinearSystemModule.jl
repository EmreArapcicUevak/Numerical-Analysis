module LinearSystemModule

    using MatrixAnalysisModule, ElemenetaryMatrixOperationsModule
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

    export LUFactorization!
    function LUFactorization!(A :: Matrix{T}) <: Matrix{T} where T <: validType
        rows, columns = size(A)
        L = GetIdentityMatrix(T, rows)

        for i ∈ 1:columns
                for j ∈ i+1:rows
                    c = a[j,i] / a[i,i]
                    AddRowMultiple!(A, j, i, -c)
                    AddRowMultiple!(L, j, i, c)
                end
        end
        
        return L
    end
end