using Documenter
import Pkg
Pkg.activate(".")
using BondGraphs

makedocs(
    sitename="BondGraphs.jl",
    modules = [BondGraphs],
    pages = [
        "Introduction" => "index.md",
        "Examples" => "examples.md",
        "One Ports" => "OnePorts.md",
        "Sources" => "Sources.md",
        "Two Ports" => "TwoPorts.md",
        "Junctions" => "Junctions.md",
        "Multi-Ports" => "MultiPorts.md",
        "Transfer Functions" => "TransferFunctions.md",
        "I/O" => "Import_Export.md",
        ]
        )
