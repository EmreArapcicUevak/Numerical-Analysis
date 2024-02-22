module BisectionMethodModule
    
export findRoot
export smartBisectionMethod    

rootExistance(func :: Function, a :: Float64, b :: Float64) = (sign ∘ func)(a) != (sign ∘ func)(b) # Bisection technique/Method {not as fast but it is reliable}
# Avoid multiplying because of overflow problems with big numbers

function findRoot(f :: Function, a :: Float64, b :: Float64)
    @assert a < b 

    if f(a) == 0
        return a
    elseif f(b) == 0
        return b
    elseif rootExistance(f,a,b) == false 
        return nothing
    end

    local itterations, c = 0, 0
    while b-a > 1e-10
        itterations += 1
        c = a + (b-a)/2

        if isapprox(f(c), 0) # isapprox is a function that checks if two numbers are approximately equal, helps with floating point errors
            return c, itterations
        elseif rootExistance(f, a, c)
            b = c
        else
            a = c
        end
    end

    return c, itterations
end

function smartBisectionMethod(f :: Function, a :: Float64, b :: Float64; maxItterations = 1000, shift = (x,y) -> (y-x) / 2)
    @assert a < b

    for i ∈ 1:maxItterations
        local result = findRoot(f,a,b)

        if result ≡ nothing
            local shiftAmount = shift(a,b)
            a, b = a - shiftAmount, b + shiftAmount
        else
            return result
        end
    end

    return nothing
end

end