import Pkg; Pkg.add("Documenter")

using Documenter
using BondGraphs
Pkg.add("ModelingToolkit")
using ModelingToolkit
Pkg.add("MetaGraphs")
using MetaGraphs

makedocs(
    sitename="BondGraphs.jl",
    modules = [BondGraphs],
    pages = [
        "Introduction" => "index.md",
        "Examples" => "examples.md",
        "Bond Graphs" => "BondGraphs.md",
        "One Ports" => "OnePorts.md",
        "Sources" => "Sources.md",
        "Two Ports" => "TwoPorts.md",
        "Junctions" => "Junctions.md",
        "Multi-Ports" => "MultiPorts.md",
        "Transfer Functions" => "TransferFunctions.md",
        "I/O" => "Import_Export.md",
        ]
    )
       
deploydocs(
    repo = "github.com/TRACER-LULab/BondGraphs.jl.git",
    devbranch = "main"
)

