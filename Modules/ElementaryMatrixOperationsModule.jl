module ElemenetaryMatrixOperationsModule
    using MatrixAnalysisModule
    export SwapRows!, AddRowMultiple!, MultiplyRow!

    validType = Union{Rational, AbstractFloat}

    function SwapRows!(A :: Matrix{T}, i :: Int, j :: Int) where T <: validType
        @assert i > 0 && i <= size(A,1) "i must be a valid row index"
        @assert j > 0 && j <= size(A,1) "j must be a valid row index"

        A[i,:], A[j,:] = A[j,:], A[i,:]
    end

    function AddRowMultiple!(A :: Matrix{T}, i :: Int, j :: Int, α :: T) where T <: validType
        @assert i > 0 && i <= size(A,1) "i must be a valid row index"
        @assert j > 0 && j <= size(A,1) "j must be a valid row index"

        A[i,:] += α * A[j,:]
    end

    function MultiplyRow!(A :: Matrix{T}, i :: Int, α :: T) where T <: validType
        @assert i > 0 && i <= size(A,1) "i must be a valid row index"

        A[i,:] *= α
    end
end