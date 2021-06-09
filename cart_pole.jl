using BondGraphs
using DifferentialEquations
using ModelingToolkit
using Plots

# Setup Empty BondGraph
@variables t x(t) v(t) θ(t) ω(t)
cart_pole = BondGraph(t)

