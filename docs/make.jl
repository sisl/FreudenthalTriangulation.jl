push!(LOAD_PATH, "../src/")

using Documenter, FreudenthalTriangulations

makedocs(
    modules = [FreudenthalTriangulations],
    format = Documenter.HTML(),
    sitename = "FreudenthalTriangulations.jl",
    pages = [
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "Table of Contents" => [
            "index.md",
            "install.md",
            "usage.md",
            "concepts.md"
           ]
    ]
)

deploydocs(
    repo = "github.com/sisl/FreudenthalTriangulations.jl.git",
)
