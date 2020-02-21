
using Pkg
Pkg.add("CSV")
Pkg.add("SpecialFunctions")

using CSV
using LinearAlgebra
using SpecialFunctions


X = CSV.read("C:/Users/dhojo/Desktop/X.csv")
y = CSV.read("C:/Users/dhojo/Desktop/y.csv")
n = size(X)[1]

X = convert(Matrix,X)
y = convert(Matrix,y)

r_0 = 0.06
s_0 = 6
Q_0 = [r_0 0;0 s_0]
a_0 = 6
b_0 = 600^2

## the exact value of marginal likelihood

M = X' * X + Q_0
R = I(n) - X * inv(M) * X'
ML = π^(-n/2) * b_0^(a_0/2) * (gamma((n+a_0)/2) / gamma(a_0/2)) *
  (det(Q_0)^(1/2) / det(M)^(1/2)) * ((y'*R*y)[1,1]+b_0)^(-(n+a_0)/2)
log(ML)

logML = -(n/2)*log(pi) + (a_0/2)*log(b_0) + log( gamma((n+a_0)/2)/gamma(a_0/2) ) +
  (1/2)*log( det(Q_0)/det(M) ) - ((n+a_0)/2)*log((y'*R*y)[1,1]+b_0)
logML






