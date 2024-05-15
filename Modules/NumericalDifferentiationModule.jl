module NumericalDifferentiationModule
    export ForwardDifference
    function ForwardDifference(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (f(x + h) - f(x)) / h
    end

    export BackwardDifference
    function BackwardDifference(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (f(x) - f(x - h)) / h
    end

    export CentralDifference
    function CentralDifference(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (f(x + h) - f(x - h)) / (2h)
    end

    export ThreePointForwardDifference
    function ThreePointForwardDifference(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (-3f(x) + 4f(x + h) - f(x + 2h)) / (2h)
    end

    export ThreePointBackwardDifference
    function ThreePointBackwardDifference(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (3f(x) - 4f(x - h) + f(x - 2h)) / (2h)
    end

    export CentralDifferenceSecondOrder
    function CentralDifferenceSecondOrder(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (f(x + h) - 2f(x) + f(x - h)) / h^2
    end

    export ThreePointForwardDifferenceSecondOrder
    function ThreePointForwardDifferenceSecondOrder(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (f(x) - 2f(x + h) + f(x + 2h)) / h^2
    end

    export ThreePointBackwardDifferenceSecondOrder
    function ThreePointBackwardDifferenceSecondOrder(f :: Function, x :: Real; h = 1e-6 :: Real) :: Real
        return (f(x - 2h) - 2f(x - h) + f(x)) / h^2
    end
end