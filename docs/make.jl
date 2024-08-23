using Documenter, PauliStrings
makedocs(
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "PauliStrings.jl",
    authors  = "Nicolas Loizeau",
    pages = [
        "Getting started" => "home.md",
        "Documentation" => "documentation.md",
        "Index" => "index.md",
    ]
)
