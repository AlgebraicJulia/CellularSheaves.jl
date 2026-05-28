module Style

using Plots

export set_default_style!, set_large_style!

function set_default_style!()
    default(
        thickness_scaling=1.7,
        legendfontsize=7,
        guidefontsize=10,
        tickfontsize=8,
    )
    return nothing
end

function set_large_style!()
    default(
        thickness_scaling=2.0,
        legendfontsize=9,
        guidefontsize=12,
        tickfontsize=10,
    )
    return nothing
end

end # module Style
