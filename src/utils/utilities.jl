# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
utilities.jl — Pure utility routines ported from CSHORE_USACE_hardbottom.f
==============================================================================#

"""
    trapz_integrate(y, dx) -> Float64

Trapezoidal integration with uniform spacing `dx`. Accumulates
(y[j] + y[j+1])*dx/2 over all intervals.
"""
function trapz_integrate(y::AbstractVector{Float64}, dx::Float64)
    n = length(y)
    n < 2 && return 0.0
    s = 0.0
    @inbounds for j in 1:n-1
        s += (y[j] + y[j+1])
    end
    return 0.5 * dx * s
end

"""
    smooth_tridiagonal!(v, nsmooth=1)

In-place 3-point weighted smoothing: `v[j] ← (v[j-1] + 2·v[j] + v[j+1])/4`.
Boundary nodes are preserved. `nsmooth` repeats the pass.

Used on `DELZB` inside the Exner solver to stabilize the bed update.
"""
function smooth_tridiagonal!(v::AbstractVector{Float64}, nsmooth::Int=1)
    n = length(v)
    n < 3 && return v
    tmp = similar(v)
    for _ in 1:nsmooth
        @inbounds begin
            tmp[1] = v[1]
            tmp[n] = v[n]
            for j in 2:n-1
                tmp[j] = 0.25 * (v[j-1] + 2 * v[j] + v[j+1])
            end
            copyto!(v, tmp)
        end
    end
    return v
end

"""
    extrapolate_boundary!(v)

Linearly extrapolates the two boundary nodes of `v` from the interior slope:
`v[1] = 2·v[2] - v[3]`, `v[n] = 2·v[n-1] - v[n-2]`.

Used after smoothing to preserve the edge values of DELZB.
"""
function extrapolate_boundary!(v::AbstractVector{Float64})
    n = length(v)
    if n ≥ 3
        v[1] = 2 * v[2] - v[3]
        v[n] = 2 * v[n-1] - v[n-2]
    end
    return v
end

"""
    interp1(x, y, xq) -> Float64

1D linear interpolation of `y(x)` at query point `xq`. Extrapolates with the
boundary values (flat extrapolation).

Assumes `x` is sorted ascending.
"""
function interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("x and y must be same length"))
    if n == 0
        return 0.0
    elseif n == 1
        return Float64(y[1])
    end
    # Clamp (flat) extrapolation
    if xq ≤ x[1]
        return Float64(y[1])
    elseif xq ≥ x[n]
        return Float64(y[n])
    end
    # Binary search for bracketing interval
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if x[mid] ≤ xq
            lo = mid
        else
            hi = mid
        end
    end
    x1, x2 = x[lo], x[hi]
    y1, y2 = y[lo], y[hi]
    t = (xq - x1) / (x2 - x1)
    return Float64(y1 + t * (y2 - y1))
end

"""
    time_series_interp(t_in, v_in, t) -> Float64

Piecewise-constant-within-interval time-series interpolation used for boundary
conditions (HRMSBC, TPBC, etc.). For `t` in `[t_in[i], t_in[i+1])`, returns
`v_in[i]` (the left endpoint value).

Out-of-range queries return the nearest endpoint value.
"""
function time_series_interp(t_in::AbstractVector{<:Real},
                            v_in::AbstractVector{<:Real},
                            t::Real)
    n = length(t_in)
    n == length(v_in) || throw(DimensionMismatch("t_in and v_in must match"))
    n == 0 && return 0.0
    if t ≤ t_in[1]
        return Float64(v_in[1])
    elseif t ≥ t_in[n]
        return Float64(v_in[n])
    end
    # Binary search for i such that t_in[i] ≤ t < t_in[i+1]
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        if t_in[mid] ≤ t
            lo = mid
        else
            hi = mid
        end
    end
    return Float64(v_in[lo])
end

"""
    erfcc(x)

Complementary error function. Thin wrapper around `SpecialFunctions.erfc`.
"""
erfcc(x::Real) = erfc(x)
