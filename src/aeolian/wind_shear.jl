# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
wind_shear.jl — Wind-flow-over-topography shear-stress perturbation.

Spatially varies the friction velocity over a non-flat bed so the aeolian
kernel sees the right local u* — bigger over the dune crest (speed-up),
smaller in the lee shadow (separation), and unchanged in the wet zone
(masked out so submerged bars don't generate spurious shear bumps).

Two pieces:

  1. **Kroy / Sauermann / Herrmann (2002) "minimal-model" perturbation**:
     an analytical Hilbert-transform of the local bed slope plus a
     local-slope term:

         τ/τ₀(x) = 1 + α [ H[∂z/∂x](x) + β · ∂z/∂x(x) ]

     The Hilbert transform captures the upstream-asymmetric pressure
     perturbation that decelerates flow on the stoss face approach,
     accelerates over the crest, and decelerates again past the brink;
     the local-slope term tunes the magnitude of stoss-vs-lee asymmetry.
     Inputs: bed elevation z(x) on a uniform-dx grid. Output: τ/τ₀(x).

     Two implementations are provided, selectable via
     `WindShearConfig.kroy_solver`:

       :direct (default) — discrete real-space convolution, O(N²):
           `H[f](x_i) ≈ (1/π) Σ_{j ≠ i} f(x_j) · dx / (x_i - x_j)`.
         Direct convolution has no boundary artifacts since it never
         assumes periodicity.

       :fft — O(N log N) via the spectral identity that the Hilbert
         transform is multiplication by `−i · sign(k)` in Fourier
         space. Uses `FFTW.jl` (real-valued FFT). For high-resolution
         decadal runs (N ≳ 500), this is materially faster. Zero-padded
         (non-periodic) FFT can introduce edge ringing in the slope
         field; we mitigate by linearly extending the input outside the
         domain before the FFT, which removes the boundary discontinuity
         that drives Gibbs artifacts.

  2. **Lee separation mask**: project a constant-slope plane downwind of
     every landward-facing brink (where the bed transitions from rising
     to falling). Cells under that plane are inside the recirculation
     bubble; their τ is collapsed to 0.

The output `τ/τ₀(x)` is multiplied by the uniform u* in the aeolian
kernel, giving a per-cell `u*(x) = u*_uniform · √(τ/τ₀)`.

Reference for the Kroy perturbation:
- Kroy, K., Sauermann, G., & Herrmann, H.J. (2002). "Minimal model for
  sand dunes." Phys. Rev. Lett., 88, 054301.
- Roelvink, D. — DUNA implementation: https://github.com/danoroelvink/duna
==============================================================================#

"""
    kroy_shear_perturbation(x::Vector, zb::Vector; alpha=3.0, beta=0.2,
                             clamp_floor=0.1) -> Vector

Compute τ/τ₀(x) at each cell from the Kroy-Sauermann-Herrmann minimal
model. `x` and `zb` must be co-located, equispaced (dx = constant), and
have the same length.

- `alpha` — overall amplitude of the perturbation (DUNA default 3.0).
  Tunes how much the shear responds to bed-slope features.
- `beta` — local-slope contribution weight (DUNA default 0.2). Controls
  stoss-vs-lee asymmetry: positive β biases the perturbation onshore
  of the crest.
- `clamp_floor` — minimum allowed τ/τ₀ (DUNA default 0.1). Prevents
  unphysical negative shear in deep lee zones; the lee separation
  mask handles the actual bubble.

Returns a vector of length `length(x)` with τ/τ₀ at each cell.
"""
function kroy_shear_perturbation(x::Vector{Float64}, zb::Vector{Float64};
                                  alpha::Float64 = 3.0,
                                  beta::Float64 = 0.2,
                                  clamp_floor::Float64 = 0.1,
                                  solver::Symbol = :direct)
    solver in (:direct, :fft) || throw(ArgumentError(
        "kroy_shear_perturbation: solver must be :direct or :fft, got :$solver"))
    n = length(x)
    n == length(zb) || throw(DimensionMismatch(
        "x and zb must have the same length (got $n vs $(length(zb)))"))
    n < 3 && return ones(Float64, n)
    dx = x[2] - x[1]
    dx > 0 || throw(ArgumentError("x must be ascending and equispaced (got dx=$dx)"))

    # Local bed slope ∂z/∂x via central differences (one-sided at boundaries)
    dzdx = Vector{Float64}(undef, n)
    @inbounds begin
        dzdx[1] = (zb[2] - zb[1]) / dx
        for i in 2:(n - 1)
            dzdx[i] = (zb[i + 1] - zb[i - 1]) / (2 * dx)
        end
        dzdx[n] = (zb[n] - zb[n - 1]) / dx
    end

    # Hilbert transform of dzdx — choice of solver
    Hdz = solver == :direct ?
        _hilbert_direct(dzdx) :
        _hilbert_fft(dzdx)

    # Kroy perturbation
    tau_ratio = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        tau_ratio[i] = max(clamp_floor, 1.0 + alpha * (Hdz[i] + beta * dzdx[i]))
    end
    return tau_ratio
end

"""
    _hilbert_direct(f) -> Vector

Discrete Hilbert transform via direct O(N²) convolution:
`H[f](x_i) ≈ (1/π) · Σ_{j ≠ i} f(x_j) / (i − j)`
(the `dx` cancels because `(x_i - x_j) = (i-j)·dx` and we have a `dx`
in the numerator from the Σ over cells).
"""
function _hilbert_direct(f::Vector{Float64})
    n = length(f)
    inv_pi = 1.0 / π
    out = zeros(Float64, n)
    @inbounds for i in 1:n
        s = 0.0
        for j in 1:n
            i == j && continue
            s += f[j] / (i - j)
        end
        out[i] = inv_pi * s
    end
    return out
end

"""
    _hilbert_fft(f) -> Vector

Discrete Hilbert transform via FFT, O(N log N). Uses the spectral
identity that the Hilbert transform corresponds to multiplication by
`−i · sign(k)` in Fourier space. Real-input rfft / irfft for efficiency.

To suppress edge ringing (Gibbs artifacts from non-periodic input), the
slope field is symmetrically reflected before the transform, then
cropped back to the original length. The reflection makes the
extended signal continuous at both boundaries — the FFT only "sees"
a smoothly periodic signal, so no spurious high-k content gets
generated by the wrap-around discontinuity.
"""
function _hilbert_fft(f::Vector{Float64})
    n = length(f)
    n < 4 && return _hilbert_direct(f)        # FFT not worth it for tiny n

    # Even-symmetric reflection: extend f → [f; reverse(f)] (length 2n).
    # This is C¹-smooth at both endpoints (since reverse(f) ends with f[1]
    # and starts with f[end], no jumps), so the periodic FFT doesn't
    # create boundary spikes.
    g = vcat(f, reverse(f))
    N = length(g)

    F = FFTW.rfft(g)
    # Multiply by -i · sign(k) — the spectral Hilbert kernel.
    # rfft frequencies: k = 0, 1, 2, ..., N/2 (all non-negative).
    # sign(0) = 0, sign(k>0) = +1.
    F[1] = 0.0 + 0.0im                         # DC: zeroed
    @inbounds for k in 2:length(F)
        F[k] = -im * F[k]
    end
    g_h = FFTW.irfft(F, N)
    return g_h[1:n]                            # crop back to original length
end

"""
    aeolis_kroy_alpha_beta(L, z0; karman=0.4) -> (α, β)

Compute Kroy 2002 self-consistent `(α, β)` from a feature length-scale
`L` (m) and roughness `z0` (m). In dimensionless form with `Φ = L/z₀`
and `φ = ℓ_inner/z₀`, the inner layer satisfies the implicit relation

    φ · log(φ) = 2κ² · Φ

(solved here by Newton-Raphson). Then

    fact = 1 + log(φ) + 2·log(π/2) + 4·γ_E
    α = (log(Φ² / log Φ))² / (2 · log(φ)³) · fact
    β = π / fact

For typical coastal foredune scales (L ≈ 10–30 m, z₀ ≈ 1e-5–1e-4 m),
`α ≈ 3–5` and `β ≈ 0.20–0.25`. The self-consistent formulation gives a
principled parameter dependence on the feature scale rather than a single
calibrated number, useful when running across very different dune
geometries with the same model.
"""
function aeolis_kroy_alpha_beta(L::Float64, z0::Float64; karman::Float64 = 0.4)
    L > 0 && z0 > 0 || throw(ArgumentError(
        "L and z0 must be > 0 (got L=$L, z0=$z0)"))
    γ_em = 0.5772156649015329
    Φ = L / z0
    log_Φ = log(Φ)
    rhs = 2 * karman * karman * Φ

    # Newton solve for φ: f(φ) = φ - rhs/log(φ) = 0
    # Initial guess from the regime where φ·log(φ) ≈ rhs:
    φ = max(rhs / log_Φ, 2.0)
    for _ in 1:80
        lp = log(φ)
        f  = φ - rhs / lp
        f′ = 1.0 + rhs / (φ * lp * lp)
        Δ  = f / f′
        φ_new = φ - Δ
        if !isfinite(φ_new) || φ_new <= 1.0
            φ_new = 0.5 * (φ + 1.001)         # safeguard
        end
        if abs(φ_new - φ) < 1e-10 * abs(φ_new)
            φ = φ_new; break
        end
        φ = φ_new
    end

    log_φ = log(φ)
    fact = 1.0 + log_φ + 2.0 * log(0.5π) + 4.0 * γ_em
    α = ( log(Φ * Φ / log_Φ) )^2 / (2 * log_φ^3) * fact
    β = π / fact
    return α, β
end

"""
    gaussian_smooth_mass_conserving(h, sigma, dx) -> Vector

Convolve `h` with a discrete Gaussian (σ in metres, ±3σ truncated,
zero-padded boundaries), then rescale so the smoothed area equals the
input area. Used to filter sub-meter Hilbert-transform noise out of the
τ/τ₀ profile while preserving the macroscopic crest speed-up. Returns
a copy; does not modify input.
"""
function gaussian_smooth_mass_conserving(h::Vector{Float64}, sigma::Float64,
                                          dx::Float64)
    n = length(h)
    sigma > 0 || return copy(h)
    half_w = max(1, ceil(Int, 3 * sigma / dx))
    inv2sig2 = 1.0 / (2 * sigma * sigma)
    kernel = Float64[exp(-((i * dx)^2) * inv2sig2) for i in -half_w:half_w]
    kernel ./= sum(kernel)
    out = zeros(Float64, n)
    @inbounds for i in 1:n
        s = 0.0
        for k in -half_w:half_w
            ji = i + k
            if 1 <= ji <= n
                s += kernel[k + half_w + 1] * h[ji]
            end
        end
        out[i] = s
    end
    area_in  = sum(h)
    area_out = sum(out)
    if area_out > 1e-15
        out .*= area_in / area_out
    end
    return out
end

"""
    lee_separation_mask(x::Vector, zb::Vector; slope=0.4,
                        min_brink_drop=0.0) -> Vector{Bool}

Identify cells inside lee-separation bubbles. For each "brink" (cell
where the bed transitions from rising landward to falling landward),
project a downwind plane `z_lee(x') = zb_brink - slope · (x' - x_brink)`
landward; cells where `zb < z_lee` are inside the bubble and have their
shear stress effectively zeroed.

- `slope` — recirculation-bubble slope (DUNA default 0.4 = ~22°). Steeper
  values produce shorter, deeper bubbles; gentler values give longer
  shadows.
- `min_brink_drop` (m) — only consider brinks where the bed drops by at
  least this much within ~5 cells downwind. Filters out tiny ripples.
"""
function lee_separation_mask(x::Vector{Float64}, zb::Vector{Float64};
                              slope::Float64 = 0.4,
                              min_brink_drop::Float64 = 0.0,
                              wind_direction::Symbol = :onshore)
    n = length(x)
    mask = falses(n)
    n < 3 && return mask
    dx = x[2] - x[1]
    look = max(1, ceil(Int, 5))     # 5-cell sniff window for "real" brinks

    # For onshore wind the brinks are landward-facing peaks: cells where
    # the bed rises seaward of and falls landward of the cell. We sweep
    # landward (j increasing) and at each candidate brink project the
    # lee plane.
    if wind_direction == :onshore
        @inbounds for i in 2:(n - 1)
            rising_seaward = zb[i] > zb[i - 1]
            falling_landward = zb[i] > zb[min(i + 1, n)]
            if !(rising_seaward && falling_landward); continue; end
            # Confirm a real drop within `look` cells
            j_end = min(n, i + look)
            if zb[i] - zb[j_end] < min_brink_drop; continue; end
            # Project the lee plane landward from the brink
            zb_brink = zb[i]
            x_brink  = x[i]
            for j in (i + 1):n
                z_lee = zb_brink - slope * (x[j] - x_brink)
                if zb[j] < z_lee
                    mask[j] = true
                else
                    break               # plane has hit the bed → bubble closed
                end
            end
        end
    else
        # Offshore wind: mirror image. Brinks are seaward-facing peaks.
        @inbounds for i in 2:(n - 1)
            rising_landward = zb[i] > zb[min(i + 1, n)]
            falling_seaward = zb[i] > zb[max(i - 1, 1)]
            if !(rising_landward && falling_seaward); continue; end
            j_end = max(1, i - look)
            if zb[i] - zb[j_end] < min_brink_drop; continue; end
            zb_brink = zb[i]
            x_brink  = x[i]
            for j in (i - 1):-1:1
                z_lee = zb_brink - slope * (x_brink - x[j])
                if zb[j] < z_lee
                    mask[j] = true
                else
                    break
                end
            end
        end
    end
    return mask
end

"""
    apply_wet_mask!(zb_eff::Vector, swl::Float64)

In-place: replace any cell where `zb_eff < swl` with `swl`. Submerged
bars and scour holes are excluded from the shear perturbation calculation
so the wind sees an effectively flat water surface at the SWL elevation
over wet cells.
"""
function apply_wet_mask!(zb_eff::Vector{Float64}, swl::Float64)
    @inbounds for i in eachindex(zb_eff)
        if zb_eff[i] < swl
            zb_eff[i] = swl
        end
    end
    return zb_eff
end

"""
    compute_wind_shear!(state, config, l, swl) -> nothing

Top-level wind-shear-over-topography update for cross-shore line `l`.
When `iwindshear == 1`, builds an effective topography (bathymetry
masked at SWL), runs the Kroy perturbation, applies the lee separation
mask, and stores per-cell `τ/τ₀` in `state.tau_perturbation` plus a
boolean `state.lee_zone` mask. Both are read by the aeolian kernel
when computing local `u*(x) = u*_uniform · √(τ/τ₀)`.

When `iwindshear == 0`, fills `tau_perturbation` with 1.0 and clears
`lee_zone`, so downstream code sees a uniform u* without needing to
special-case anything.
"""
function compute_wind_shear!(state::CshoreState, config::CshoreConfig,
                              l::Int, swl::Float64;
                              wind_direction::Symbol = :onshore)
    jmax = state.jmax[l]

    # Default: no perturbation
    @inbounds for j in 1:jmax
        state.tau_perturbation[j] = 1.0
        state.lee_zone[j] = false
    end

    config.options.iwindshear == 1 || return nothing
    cfg_ws = config.windshear
    cfg_ws === nothing && return nothing

    # Build the effective topography to feed the Kroy solver
    x = view(state.xb, 1:jmax)
    zb_eff = Vector{Float64}(undef, jmax)
    @inbounds for j in 1:jmax
        zb_eff[j] = state.zb[j, l]
    end
    if cfg_ws.mask_bathymetry
        apply_wet_mask!(zb_eff, swl)
    end

    # Choose α, β based on shear_method.
    α, β = if cfg_ws.shear_method == :aeolis_kroy
        # Self-consistent (α, β) from Kroy 2002 implicit solve. Roughness
        # z0 is taken from the surface-D50 weighting (the same path the
        # aeolian kernel uses), with a fallback to 2/30 · sediment.d50.
        nf = nfractions(config.multifraction)
        d50_local = if nf > 0
            d_k_sum = 0.0
            w_sum   = 0.0
            for k in 1:nf
                w = sum(view(state.active_frac, 1:jmax, k)) / max(jmax, 1)
                d_k_sum += w * config.multifraction.grain_sizes[k]
                w_sum   += w
            end
            w_sum > 0 ? d_k_sum / w_sum : config.sediment.d50
        else
            config.sediment.d50
        end
        z0 = max(2.0 / 30.0 * d50_local, 1e-6)
        α_ae, β_ae = aeolis_kroy_alpha_beta(cfg_ws.aeolis_length_scale, z0)
        α_ae, β_ae
    else
        # User-specified (α, β).
        cfg_ws.kroy_alpha, cfg_ws.kroy_beta
    end

    # Kroy shear perturbation
    tau = kroy_shear_perturbation(collect(x), zb_eff;
                                   alpha = α, beta = β,
                                   clamp_floor = cfg_ws.tau_clamp_floor,
                                   solver = cfg_ws.kroy_solver)

    # Optional Gaussian smoothing of τ/τ₀ (filters Hilbert-transform
    # grid-noise spikes; mass-conserving rescale).
    if cfg_ws.tau_smooth_sigma > 0
        dx = state.xb[2] - state.xb[1]
        tau = gaussian_smooth_mass_conserving(tau, cfg_ws.tau_smooth_sigma, dx)
        # The smoothing can drift below the floor at deep lee zones; reapply.
        @inbounds for j in 1:length(tau)
            tau[j] = max(tau[j], cfg_ws.tau_clamp_floor)
        end
    end

    @inbounds for j in 1:jmax
        state.tau_perturbation[j] = tau[j]
    end

    # Lee separation mask — uses the *real* bathymetry, not the wet-masked
    # version, since the lee bubble is a feature of the actual brink shape.
    zb_real = Vector{Float64}(undef, jmax)
    @inbounds for j in 1:jmax
        zb_real[j] = state.zb[j, l]
    end
    mask = lee_separation_mask(collect(x), zb_real;
                                slope = cfg_ws.lee_slope,
                                min_brink_drop = cfg_ws.min_brink_drop,
                                wind_direction = wind_direction)
    @inbounds for j in 1:jmax
        state.lee_zone[j] = mask[j]
        if mask[j]
            state.tau_perturbation[j] = 0.0
        end
    end
    return nothing
end
