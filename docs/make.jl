using SystemBondGraphs
using Documenter
# using DocumenterCitations

# bib = CitationBibliography(joinpath(@__DIR__, "src", "paper.bib"); style = :numeric)

DocMeta.setdocmeta!(SystemBondGraphs, :DocTestSetup, :(using SystemBondGraphs); recursive = true)

makedocs(;
    # plugins = [bib],
    modules = [SystemBondGraphs],
    authors = "Carson Farmer <59753859+cfarm6@users.noreply.github.com> and contributors",
    repo = "https://github.com/TRACER-LULab/SystemBondGraphs.jl/blob/{commit}{path}#{line}",
    sitename = "SystemBondGraphs.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://TRACER-LULab.github.io/SystemBondGraphs.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", "API" => "API.md", "Example" => "example.md"],
)

deploydocs(; repo = "github.com/TRACER-LULab/SystemBondGraphs.jl", devbranch = "main")
