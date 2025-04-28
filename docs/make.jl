using Documenter, PauliStrings
using PauliStrings.Circuits




readme_str = read(joinpath(@__DIR__, "..", "README.md"), String)
write(
    joinpath(@__DIR__, "src", "index.md"),
    replace(readme_str, "./docs/src/assets/" => "assets/"),
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
            "Circuits" => "circuits.md"],
        "Docstrings" => "docstrings.md",
        "Index" => "docstrings_index.md"]
)

deploydocs(
    repo="github.com/nicolasloizeau/PauliStrings.jl.git",
    devbranch="main",
    branch="gh-pages",
)
