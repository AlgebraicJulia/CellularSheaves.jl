############################################################################################
# accum_body
############################################################################################

function accum_body(nout, nred, cbase::Symbol, bbase::Symbol, aidx, ta)
    Δ = [Symbol(:Δ, k) for k in 1:nout]

    init  = [:($(Δ[k]) = zero(promote_eltype(A, B))) for k in 1:nout]
    accum = [:($(Δ[k]) = muladd(wrapadj(A[bstrt + $(aidx(k))], $ta), Bx, $(Δ[k]))) for k in 1:nout]
    store = [:(C[$cbase + $(k - 1)] = muladd(α, $(Δ[k]), C[$cbase + $(k - 1)])) for k in 1:nout]

    return quote
        $(init...)

        @inbounds for t in 1:$nred
            Bx = B[$bbase + t - 1]

            $(accum...)
        end

        @inbounds begin
            $(store...)
        end
    end
end

include("gemv.jl")
include("gemm.jl")
