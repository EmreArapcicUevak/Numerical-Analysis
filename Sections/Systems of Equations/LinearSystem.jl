push!(LOAD_PATH, pwd()*"/Modules/")
using LinearSystemModule
using BisectionMethodModule, Revise, Plots, ShowPointModule, HornerMethodModule, NewtonMethodModule, LinearSystemModule
plotlyjs() # Set PlotlyJS as the backend for Plots.jl

A = [1 4/3 5/3; 0 1 11/10; 0 0 1]
b = [10/3, 13/10, 3]

SolveUpperDiagonal(A, b)
SolveUpperDiagonal([1 4//3 5//3; 0 1 4//10; 0 0 1], [10//3, 13//10, 3])

A = [-2. 0 0; 4. -1 0; 0 4 -1]
b = Float64[2, 1, 1]
SolveLowerDiagonal(A, b)

A = Float64[-1 0 3; 2 1 3; 1 1 2]
result = PLUFactorization(A)
A
result.U
result.L
result.P

A = Rational[3 -7 -2 2; -3 5 1 0; 6 -4 0 -5; -9 5 -5 12]
result = PLUFactorization(A)
A
result.U
result.L
result.P


A = Rational[2 -4 -2 3; 6 -9 -5 8; 2 -7 -3 9; 4 -2 -2 -1; -6 3 3 4]
result = PLUFactorization(A)
A
result.U
result.L
result.L * result.U

A = Rational[2 3 1; 4 1 4; 3 4 6]
result = PLUFactorization(A)
A
result.U
result.L
result.P

b = Rational[-4, 9, 0]
SolvePLU(result.P, result.L, result.U, b)
y = SolveLowerDiagonal(result.L, b)
x = SolveUpperDiagonal(result.U, y)

A * x 
result.L * result.U

A = Float64[1 1 1; 1 2 2; 1 1 2]
result = PLUFactorization(A)
result.U
result.L 
result.L * result.U

A = Rational[3 -7 -2 2; -3 5 1 0; 6 -4 0 -5; -9 5 -5 12]
result = PLUFactorization(A)
A
result.U
result.L
result.P

b = Rational[-9, 5, 7, 11]
x = SolvePLU(result.P, result.L, result.U, b)
result.L * result.U * x 

A = Float64[-1 0 3; 2 1 3; 1 1 2]
result = PLUFactorization(A)
result.U
result.L
result.P

A = Rational[1 1 -1; 1 -2 3; 2 3 1]
b = Rational[4, -6, 7]
result = PLUFactorization(A)
result.U
result.L
result.P
typeof(result)
x = SolvePLU(result, b)

A = Rational[10 -7 0; -3 2 6; 5 -1 5]
b = Rational[1,2,3]
result = PLUFactorization(A)
x = SolvePLU(result, b)

A = Rational[2 1 -1; -3 -1 2; -2 1 2]
b = Rational[1,2,3]
result = PLUFactorization(A)
result.U
result.L
result.P
x = SolvePLU(result, b)
x

A = Float64[21 -10 -1 0; -10 21 -10 -1; -1 -10 21 -10; 0 -1 -10 21]
b = Float64[50, 0, 0, 0]
result = PLUFactorization(A)
x = SolvePLU(result, b)