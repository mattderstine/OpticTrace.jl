import Pkg
Pkg.develop(Pkg.PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using OpticTrace
using Documenter

makedocs(;
    modules=[OpticTrace],
    authors="Matthew Derstine",
    sitename="OpticTrace.jl",
    format=Documenter.HTML(; prettyurls=false, size_threshold_ignore=["api.md"]),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)
