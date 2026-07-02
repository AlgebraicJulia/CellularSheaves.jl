const BRAILLE_BLOCKS = (0x2801, 0x2802, 0x2804, 0x2840, 0x2808, 0x2810, 0x2820, 0x2880)

const AbstractScalar{T} = AbstractArray{T, 0}

function axpy(a, x, y)
    return muladd(a, x, y)
end

function axpy(a::Real, x::Real, y::Real)
    return fma(a, x, y)
end

function ispositive(i::I) where {I}
    return i > zero(I)
end

function two(::Type{I}) where {I}
    return one(I) + one(I)
end

function two(::I) where {I}
    return two(I)
end

function intriangle(i, j, uplo::Symbol)
    if uplo == :U
        return i <= j
    else
        return i >= j
    end
end

function intriangle(i, j, ::Val{UL}) where {UL}
    return intriangle(i, j, UL)
end

function unwrapadj(A::Transpose)
    return parent(A), Val(:T)
end

function unwrapadj(A::Adjoint)
    return parent(A), Val(:C)
end

function unwrapadj(A)
    return A, Val(:N)
end

function unwrapsym(A::Symmetric)
    return parent(A), Val(:N)
end

function unwrapsym(A::Hermitian)
    return parent(A), Val(:C)
end

function wrapsym(A, ::Val{:N}, uplo::Symbol)
    return Symmetric(A, uplo)
end

function wrapsym(A, ::Val{:C}, uplo::Symbol)
    return Hermitian(A, uplo)
end

function wrapadj(A, ::Val{TA}) where {TA}
    if TA === :N
        B =           A
    elseif TA === :T
        B = transpose(A)
    else
        B =   adjoint(A)
    end

    return B
end

function binarysearchlast(A::AbstractVector{I}, v::Integer, strt::Integer, stop::Integer) where {I <: Integer}
    return binarysearchlast(A, convert(I, v)::I, convert(I, strt)::I, convert(I, stop)::I)
end

function binarysearchlast(A::AbstractVector{I}, v::I, strt::I, stop::I) where {I <: Integer}
    lo = strt - one(I)
    hi = stop + one(I)

    while lo < hi - one(I)
        md = (lo + hi) >> 1

        if A[md] <= v
            lo = md
        else
            hi = md
        end
    end

    return lo
end

function braille_grid(io::IO, nrow, ncol)
    maxheight, maxwidth = displaysize(io)
    maxheight -= 4
    maxwidth ÷= 2

    if get(io, :limit, true)
        scaleheight = max(8, min(nrow, 2maxwidth, 4maxheight))
        scalewidth = max(4, min(ncol, 2maxwidth, 4maxheight))
    else
        scaleheight = max(8, nrow)
        scalewidth = max(4, ncol)
    end

    rowscale = max(1, scaleheight - 1) / max(1, nrow - 1)
    colscale = max(1, scalewidth - 1) / max(1, ncol - 1)

    grow = (scalewidth - 1) ÷ 2 + 4
    gcol = (scaleheight - 1) ÷ 4 + 1

    grid = FMatrix{UInt16}(undef, grow, gcol)
    grid                 .= '⠀'
    grid[1,           :] .= '⎢'
    grid[grow - 1,    :] .= '⎥'
    grid[1,           1]  = '⎡'
    grid[1,        gcol]  = '⎣'
    grid[grow - 1,    1]  = '⎤'
    grid[grow - 1, gcol]  = '⎦'
    grid[grow,        :] .= '\n'

    return grid, rowscale, colscale
end

function setbraille!(braillegrid, row, col, rowscale, colscale)
    si = round(Int, (row - 1) * rowscale + 1)
    sj = round(Int, (col - 1) * colscale + 1)

    i =  (sj - 1) ÷ 2  + 2
    j =  (si - 1) ÷ 4  + 1
    b = ((sj - 1) % 2) * 4 + ((si - 1) % 4 + 1)

    braillegrid[i, j] |= BRAILLE_BLOCKS[b]
end
