using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
# using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Symbolics
## 
@variables x f(x) q h m
lim h→0 ((-1)^q/h^q ∑_m=0^∞ (-1)^m*binomial(q, m)f(x+m*h))
