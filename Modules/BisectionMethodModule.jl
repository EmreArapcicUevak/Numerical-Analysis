module BisectionMethodModule
    
export findRoot
export smartBisectionMethod    

rootExistance(leftBound :: Real, rightBound :: Real) = sign(leftBound) ≠ sign(rightBound)

function findRoot(f :: Function, a :: Real, b :: Real; precision = 1e-15 :: Float64)
    @assert a < b 
    local f_a, f_b = f(a), f(b)

    if abs(f_a) ≤ precision
        return a
    elseif abs(f_b) ≤ precision
        return b
    elseif rootExistance(f_a, f_b) == false 
        return nothing
    end

    local itterations, c = 0, a + 0.5*(b-a)
    while b-a > precision
        itterations += 1

        f_c = f(c)
        if abs(f_c) ≤ precision # isapprox is a function that checks if two numbers are approximately equal, helps with floating point errors
            return (c = c, itterations = itterations)
        elseif rootExistance(f_a, f_c)
            b = c
        else
            a = c
        end
        c = a + 0.5*(b-a)
    end

    return (c = c, itterations = itterations)
end

function smartBisectionMethod(f :: Function, a :: Real, b :: Real, k = 2 :: Real; maxItterations = 1000 :: Integer, domain = (-Inf, Inf) :: Tuple{Real, Real}, precision = 1e-15 :: Real)
    @assert a < b
    L = b - a
    center = a + 0.5*L

    a, b = domain[1] ≤ a ≤ domain[2] ? a : domain[1], domain[1] ≤ b ≤ domain[2] ? b : domain[2]
    
    for i ∈ 1:maxItterations
        println("(a, b) = ($a, $b)")
        if rootExistance(f(a), f(b))
            println("(a, b) = ($a, $b)\nItterations taken: $i")
            return findRoot(f, a, b; precision = precision)
        else
            L *= k
            a = center - L ≥ domain[1] ? center - L : domain[1]
            b = center + L ≤ domain[2] ? center + L : domain[2]
        end
    end

    return nothing
end

end