using TSCSMethods
using Documenter

DocMeta.setdocmeta!(TSCSMethods, :DocTestSetup, :(using TSCSMethods); recursive=true)

makedocs(;
    modules=[TSCSMethods],
    authors="Eric Martin Feltham <eric.feltham@yale.edu> and contributors",
    repo=Documenter.Remotes.GitHub("emfeltham", "TSCSMethods.jl"),
    sitename="TSCSMethods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://emfeltham.github.io/TSCSMethods.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Methodology" => "methodology.md",  
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/emfeltham/TSCSMethods.jl",
)
