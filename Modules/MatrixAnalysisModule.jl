module MatrixAnalysisModule

    export IsSquare, IsDiagonal, IsUpperTriangular, IsLowerTriangular
    const validType = Union{Rational, AbstractFloat}

    function AproxZero(x :: validType)
        return abs(x) < 1e-10
    end

    function IsSquare(A :: Array{T,2}) where T <: validType
        n, m = size(A)
        return n == m
    end

    function IsDiagonal(A :: Array{T,2}) where T <: validType
        n, m = size(A)
        if n != m
            return false
        else
            for i ∈ 1:n
                for j ∈ 1:m
                    if i != j && !AproxZero(A[i,j])
                        return false
                    elseif i == j && AproxZero(A[i,j])
                        return false
                    end
                end
            end
            return true
        end
    end

    function IsUpperTriangular(A :: Array{T,2}) where T <: validType
        n, m = size(A)
        if n != m
            return false
        else
            for i ∈ 1:n
                for j ∈ 1:m
                    if i > j && !AproxZero(A[i,j])
                        return false
                    end
                end
            end

            return true
        end
    end

    function IsLowerTriangular(A :: Array{T,2}) where T <: validType
        n, m = size(A)
        if n != m
            return false
        else
            for i ∈ 1:n
                for j ∈ 1:m
                    if i < j && !AproxZero(A[i,j])
                        return false
                    end
                end
            end

            return true
        end
    end

end
