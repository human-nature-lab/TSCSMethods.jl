using TSCSMethods
using Documenter

DocMeta.setdocmeta!(TSCSMethods, :DocTestSetup, :(using TSCSMethods); recursive=true)

makedocs(;
    modules=[TSCSMethods],
    authors="Eric Martin Feltham <eric.feltham@yale.edu> and contributors",
    repo="https://github.com/human-nature-lab/TSCSMethods.jl/blob/{commit}{path}#{line}",
    sitename="TSCSMethods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://human-nature-lab.github.io/TSCSMethods.jl",
        edit_link="main",
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
    repo="github.com/human-nature-lab/TSCSMethods.jl",
)
