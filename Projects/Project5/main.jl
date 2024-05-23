push!(LOAD_PATH, pwd()*"/Modules/")
#mlθ'' = -mg sin(θ) - blθ'
# θ'' = -g/l sin(θ) - bθ'/m
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

#using Revise, Plots # Always leave this on!
#Plots.plotlyjs() # Set PlotlyJS as the backend for Plots.jl
using Javis

#### Other modules that you might need go here ####
using DifferentialEquationsModule
###################################################

#### Your code goes here ####
g = 9.81
m = 1
b = 0
l = 10
f(t, z) = [z[2], -g/l * sin(z[1]) - b/m * z[2]]

θ₀ = pi/2
θ₁ = 0

z₀ = [θ₀, θ₁]

N = 1200
T = 40

boundaries = (0.0, T)
solutions = MidpointMethod(f, z₀, N-1, boundaries)

θSolution, θPrimeSolution = Real[], Real[]
for i ∈ solutions
    push!(θSolution, i[1])
    push!(θPrimeSolution, i[2])
end

points = LinRange(boundaries[1], boundaries[2], N)
x(t) = θ₀ * cos(sqrt(g/l) * t)

#plot(points, x.(points), title = "Pendulum", label = "Exact Solution", xlabel = "t", ylabel = "θ(t)", lw = 3)
plot(points, θSolution, label = "θ Numerical Solution", lw = 2, linestyle = :dash)
plot!(points, θPrimeSolution, label = "θ' Numerical Solution", lw = 2, linestyle = :dash)



frames = length(θSolution)
scale = 0.5
video = Video(1920 * scale,1080 * scale)

function ground(args...)
    background("white")
end

trail = Point[]
function draw_pendulum_θ(video, object, frame)
    θᵢ = θSolution[frame]
    lₘ = l*20
    x,y = lₘ*sin(θᵢ), lₘ*cos(θᵢ)
    endPoint = Point(x,y)
    push!(trail, endPoint)


    sethue("black")
    line(O, endPoint, :stroke)
    circle(endPoint, m*5, :fill)
    sethue("red")
    circle.(trail, 1, :fill)
end

Background(1:frames, ground)
Object(draw_pendulum_θ)


render(video, pathname = "SwingingBall.gif")