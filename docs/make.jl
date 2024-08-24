using Documenter, PauliStrings
makedocs(
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "PauliStrings.jl",
    authors  = "Nicolas Loizeau",
    pages = [
        "Getting started" => "index.md",
        "Documentation" => "documentation.md",
        "Index" => "docstrings.md",
    ]
)

deploydocs(
    repo = "github.com/nicolasloizeau/PauliStrings.jl.git",
)
