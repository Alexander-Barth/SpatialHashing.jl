using SpatialHashing
using Documenter

DocMeta.setdocmeta!(SpatialHashing, :DocTestSetup, :(using SpatialHashing); recursive=true)

makedocs(;
    modules=[SpatialHashing],
    authors="Alexander Barth <barth.alexander@gmail.com> and contributors",
    repo="https://github.com/Alexander-Barth/SpatialHashing.jl/blob/{commit}{path}#{line}",
    sitename="SpatialHashing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Alexander-Barth.github.io/SpatialHashing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Alexander-Barth/SpatialHashing.jl",
    devbranch="main",
)
