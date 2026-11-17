using ComodoIGA
using Documenter

DocMeta.setdocmeta!(ComodoIGA, :DocTestSetup, :(using ComodoIGA); recursive=true)

makedocs(;
    modules=[ComodoIGA],
    authors="Aminofa70 <alibakhshiamin@yahoo.com> and contributors",
    sitename="ComodoIGA.jl",
    format=Documenter.HTML(;
        canonical="https://COMODO-research.github.io/ComodoIGA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/COMODO-research/ComodoIGA.jl",
    devbranch="main",
)
