using PHysicalTDA
using Documenter

DocMeta.setdocmeta!(PHysicalTDA, :DocTestSetup, :(using PHysicalTDA); recursive=true)

makedocs(;
    modules=[PHysicalTDA],
    authors="Isaac Ownby, Collin Kovacs",
    sitename="PHysicalTDA.jl",
    format=Documenter.HTML(;
    	prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://alt-f4-dev.github.io/PHysicalTDA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "tutorials.md",
        "References" => "refs.md"
    ],
)

deploydocs(;
    repo="github.com/alt-f4-dev/PHysicalTDA.jl.git",
    devbranch="main",
)
