module NumericalIntegrationModule
    export LeftRectangleRule
    function LeftRectangleRule(f :: Function, a :: Real, b :: Real; n :: Int)
        @assert n > 0 "N has to be greater than 0"
        @assert a < b "a has to be less than b"

        h = (b - a) / n
        x = LinRange(a,b, n + 1)
        y = f.(x)
        return h * sum(y[1:end-1])
    end

    export RightRectangleRule
    function RightRectangleRule(f :: Function, a :: Real, b :: Real; n :: Int)
        @assert n > 0 "N has to be greater than 0"
        @assert a < b "a has to be less than b"

        h = (b - a) / n
        x = LinRange(a,b, n + 1)
        y = f.(x)
        return h * sum(y[2:end])
    end

    export MidpointRule
    function MidpointRule(f :: Function, a :: Real, b :: Real; n :: Int)
        @assert n > 0 "N has to be greater than 0"
        @assert a < b "a has to be less than b"

        h = (b - a) / n
        x = LinRange(a,b, n + 1)
        return h * sum([f(x + h / 2) for x ∈ x[1:end-1]])
    end

    export TrapezoidalRule
    function TrapezoidalRule(f :: Function, a :: Real, b :: Real; n :: Int)
        @assert n > 0 "N has to be greater than 0"
        @assert a ≤ b "a has to be less than b"

        h = (b - a) / n
        x = LinRange(a,b, n + 1)
        y = f.(x)

        return h/2 * (y[begin] + 2*sum(y[2:end-1]) + y[end])
    end

    export SimpsonRule
    function SimpsonRule(f :: Function, a :: Real, b :: Real; n :: Int)
        @assert n > 0 "N has to be greater than 0"
        @assert a ≤ b "a has to be less than b"
        @assert n % 2 == 0 "N has to be an even number"

        h = (b-a)/n
        allX = LinRange(a,b, n+1)
        allY = f.(allX)

        return h/3 * (allY[1] + allY[end] + 4 * sum(allY[2:2:n]) + 2 * sum(allY[3:2:n-2]))
    end

    export SimpsonRule3_8
    function SimpsonRule3_8(f :: Function, a :: Real, b :: Real; n :: Int)
        @assert n > 0 "N has to be greater than 0"
        @assert a ≤ b "a has to be less than b"
        @assert n % 3 == 0 "N has to be a multiple of 3"


        h = (b-a)/n
        X = LinRange(a,b, n+1)
        Y = f.(X)

        return 3/8 * h * (Y[begin] + 2 * sum(Y[4:3:end-1]) + 3 * sum(Y[2:3:end-2]) + 3 * sum(Y[3:3:end-1]) + Y[end])
    end
end