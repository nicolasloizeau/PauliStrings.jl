using Symbolics
using Documenter, PauliStrings
using PauliStrings.Circuits
using SparseArrays
ENV["GKSwstype"] = "100"
using Plots
gr()
using LinearAlgebra

readme_str = read(joinpath(@__DIR__, "..", "README.md"), String)
index_str = replace(readme_str, "./docs/src/assets/" => "assets/")
index_str = replace(index_str, "# PauliStrings.jl" => "# Getting started")
index_str = join(split(index_str, '\n')[5:end], '\n')
write(
    joinpath(@__DIR__, "src", "index.md"),
    index_str
)

makedocs(
    format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true",
        assets=[
            "assets/custom.css",
            "assets/favicon.ico",
        ],
        analytics="G-L0F0JFQ15V"
    ),
    sitename="PauliStrings.jl",
    authors="Nicolas Loizeau",
    pages=[
        "Getting started" => "index.md",
        "Tutorials" => ["Constructing operators" => "constructing.md",
            "Lanczos" => "lanczos.md",
            "Time evolution" => "evolution.md",
            "Translation symmetry" => "translation.md",
            "Circuits" => "circuits.md",
            "Symbolics" => "symbolics.md",
            "Manipulating single strings" => "manipulating_strings.md",
            "LIOMs" => "lioms.md",
        ],
        "Docstrings" => "docstrings.md",
        "Index" => "docstrings_index.md"]
)

deploydocs(
    repo="github.com/nicolasloizeau/PauliStrings.jl.git",
    devbranch="main",
    branch="gh-pages",
)
