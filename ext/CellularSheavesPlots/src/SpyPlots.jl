module SpyPlots

using Plots
using SparseArrays

export spyplot, spyplot_pair

function spyplot(A; title::AbstractString="")
    return spy(A; title=title, legend=false)
end

function spyplot_pair(A, B; titles::Tuple{<:AbstractString,<:AbstractString}=("A", "B"))
    p1 = spy(A; title=titles[1], legend=false)
    p2 = spy(B; title=titles[2], legend=false)
    return plot(p1, p2; layout=(1, 2), size=(1000, 420))
end

end # module SpyPlots
