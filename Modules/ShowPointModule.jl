module ShowPointModule
    using Plots

    export ShowPoint

    function ShowPoint(f :: Function, x :: Real, y = 0 :: Real; label = "Point", numOfPoints = 9000, xlabel = "x", ylabel = "y" , color = :blue, title = "Plot" :: String, domain :: Tuple{Real, Real})
        @assert domain[1] ≤ x ≤ domain[2]
        
        xPoints = LinRange(domain[1], domain[2], numOfPoints)
        yPoints = f.(xPoints)

        plot(xPoints, yPoints, label = "Function", xlabel = xlabel, ylabel = ylabel, title = title, lw = 2, color = :red)
        plot!(xPoints, zeros(numOfPoints), label = "y = 0", lw = .5, color = :blue)
        plot!(zeros(numOfPoints), LinRange(minimum(yPoints), maximum(yPoints), numOfPoints), label = "x = 0", lw = .5, color = :blue)
        scatter!([x], [y], label = "Root", color = :purple)
    end

    ShowPoint(f :: Function, Point :: Tuple{Real, Real}; label = "Point", numOfPoints = 9000, xlabel = "x", ylabel = "y" , color = :blue, title = "Plot" :: String, domain :: Tuple{Real, Real}) = ShowPoint(f, Point[1], Point[2]; label = label, numOfPoints = numOfPoints, xlabel = xlabel, ylabel = ylabel, color = color, title = title, domain = domain)
end