using PIMDtools
using Documenter

DocMeta.setdocmeta!(PIMDtools, :DocTestSetup, :(using PIMDtools); recursive=true)

makedocs(;
    modules=[PIMDtools],
    authors="cometscome <cometscome@gmail.com> and contributors",
    repo="https://github.com/cometscome/PIMDtools.jl/blob/{commit}{path}#{line}",
    sitename="PIMDtools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cometscome.github.io/PIMDtools.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cometscome/PIMDtools.jl",
    devbranch="main",
)
