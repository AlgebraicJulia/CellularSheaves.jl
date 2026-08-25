# Let g be a twice-differentiable concave function
#
#   g: [0, hi] → [-∞, ∞),  g(0) ≥ 0,
#
# and let 0 ≤ τ* ≤ hi be the largest number such that
#
#   g(τ*) ≥ 0.
#
# Compute a lower bound τ ≤ τ* satisfying
#
#   τ* - τ ≤ tol (1 + hi).
#
# The argument `jet` returns a quadruple
#
#   jet(τ) = (g(τ), g'(τ), g''(τ), n(τ), c)
#
# where n(τ) is a noise bound at τ and c
# is a root candidate.
function nflast(jet, hi::T; tol::T = eps(T), maxit::Int = 53) where {T <: AbstractFloat}
    lo = zero(T); nprobe = 0

    glo, gplo, _, nlo, c = jet(lo)

    if glo >= nlo && glo >= zero(T)
        ghi, gphi, gpphi, nhi, _ = jet(hi)

        if ghi >= nhi
            lo = hi
        else
            wlo = glo

            if isfinite(nhi) && isfinite(gphi) && !iszero(gphi)
                whi = nhi / abs(gphi)
            else
                whi = zero(T)
            end

            for _ in 1:maxit
                wid = hi - lo
                #
                # lower candidate: the endgame step when hi sits inside its
                # own noise band; a candidate, then the probes, then the
                # midpoint when the face is starved (g(hi) = -Inf); otherwise
                # the chord, with chord-snap and the rescue criterion
                #
                if isfinite(nhi) && isfinite(ghi) && ghi >= -nhi
                    lon = hi - 3 * max(whi, eps(T) * (one(T) + hi)) / 2
                elseif ghi == T(-Inf)
                    #
                    # a candidate must beat the midpoint
                    #
                    if lo + (hi - lo) / 2 < c < hi - 4 * eps(T) * (one(T) + hi)
                        lon = c
                    else
                        lon = T(NaN)
                    end

                    if !(lo < lon < hi)
                        if nprobe == 0
                            lon = hi - 4 * eps(T) * (one(T) + hi)
                        elseif nprobe == 1
                            lon = hi - tol * (one(T) + hi) / 2
                        end

                        if lo < lon < hi
                            nprobe += 1
                        end
                    end

                    if !(lo < lon < hi)
                        lon = lo + (hi - lo) / 2
                    end
                else
                    lon = lo + wlo * (hi - lo) / (wlo - ghi)
                    #
                    # chord-snap: a chord target rounding onto hi certifies
                    # τ* within a float of hi — take the endgame step
                    #
                    if lon >= hi
                        lon = hi - 3 * max(whi, eps(T) * (one(T) + hi)) / 2
                    end
                    #
                    # rescue: consult the candidate only against a chord that
                    # is creeping or provably about to undershoot
                    #
                    if !isnan(c) &&
                       ((lon - lo) * 2 < hi - lo ||
                        gphi * (hi - lo) < 2 * (ghi - glo))
                        if (lo < c < hi - 4 * eps(T) * (one(T) + hi)) && c >= lon
                            lon = c
                        end
                    end
                end

                if !(lo < lon < hi)
                    lon = lo + (hi - lo) / 2
                end

                g, gp, gpp, n, k = jet(lon)

                if g >= n
                    lo, glo, gplo, wlo, c = lon, g, gp, g, k
                else
                    if isfinite(n) && isfinite(gp) && !iszero(gp)
                        whi = n / abs(gp)
                    else
                        whi = zero(T)
                    end

                    hi, ghi, gphi, gpphi, nhi = lon, g, gp, gpp, n
                    wlo /= 2
                end

                if hi - lo < max(whi, tol * (one(T) + hi))
                    break
                end
                #
                # upper candidate: single-division Halley, falling back to
                # Newton, floored by the Brent minimum step, capped by the
                # tangent from below, then the candidate, then the midpoint
                #
                hin = hi - 2 * ghi * gphi /
                           (2 * gphi * gphi - ghi * gpphi)

                if !(lo < hin < hi)
                    hin = hi - ghi / gphi
                end

                minstep = tol * (one(T) + hi) / 4

                if hi - minstep < hin < hi
                    hin = hi - minstep
                end

                if !(lo < hin < hi)
                    hin = hi
                end

                tlo = lo - glo / gplo

                if lo < tlo < hin
                    hin = tlo
                end

                if !(lo < hin < hi)
                    if lo < c < hi - 4 * eps(T) * (one(T) + hi)
                        hin = c
                    else
                        hin = T(NaN)
                    end
                end

                if !(lo < hin < hi)
                    hin = lo + (hi - lo) / 2
                end

                g, gp, gpp, n, k = jet(hin)

                if g >= n
                    lo, glo, gplo, wlo, c = hin, g, gp, g, k
                else
                    if isfinite(n) && isfinite(gp) && !iszero(gp)
                        whi = n / abs(gp)
                    else
                        whi = zero(T)
                    end

                    hi, ghi, gphi, gpphi, nhi = hin, g, gp, gpp, n
                end

                if hi - lo < max(whi, tol * (one(T) + hi))
                    break
                end
                #
                # guaranteed contraction: one bit per iteration
                #
                if 2 * (hi - lo) > wid
                    mid = lo + (hi - lo) / 2

                    g, gp, gpp, n, k = jet(mid)

                    if g >= n
                        lo, glo, gplo, wlo, c = mid, g, gp, g, k
                    else
                        if isfinite(n) && isfinite(gp) && !iszero(gp)
                            whi = n / abs(gp)
                        else
                            whi = zero(T)
                        end

                        hi, ghi, gphi, gpphi, nhi = mid, g, gp, gpp, n
                    end

                    if hi - lo < max(whi, tol * (one(T) + hi))
                        break
                    end
                end
            end
        end
    end

    return lo
end
