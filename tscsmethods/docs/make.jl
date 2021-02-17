using tscsmethods
using Documenter

DocMeta.setdocmeta!(tscsmethods, :DocTestSetup, :(using tscsmethods); recursive=true)

makedocs(;
    modules=[tscsmethods],
    authors="Eric Martin Feltham <eric.feltham@yale.edu> and contributors",
    repo="https://github.com/emfeltham/tscsmethods.jl/blob/{commit}{path}#{line}",
    sitename="tscsmethods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://emfeltham.github.io/tscsmethods.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/emfeltham/tscsmethods.jl",
)
