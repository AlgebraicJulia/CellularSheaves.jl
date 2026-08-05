using Test

@testset "docs: every @autodocs module is declared in makedocs(modules=...)" begin
    make_src = read(joinpath(@__DIR__, "..", "docs", "make.jl"), String)
    api_dir  = joinpath(@__DIR__, "..", "docs", "src", "api")
    autodoc_mods = Set{String}()
    for f in readdir(api_dir; join=true)
        endswith(f, ".md") || continue
        for m in eachmatch(r"Modules\s*=\s*\[([^\]]+)\]", read(f, String))
            for mod in split(m.captures[1], ",")
                push!(autodoc_mods, strip(mod))
            end
        end
    end
    for mod in autodoc_mods
        @test occursin(mod, make_src)
    end
end
