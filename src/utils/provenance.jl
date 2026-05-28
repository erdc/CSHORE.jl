# Developed by N. Cohn, US Army Engineer Research and Development Center —
# Coastal and Hydraulics Laboratory (2026).

#==============================================================================
provenance.jl — Sediment provenance tracking for CSHORE.jl.

Tracks which spatial source region each parcel of sediment in the active
(surface) layer originated from, throughout the simulation.

## Design

ProvenanceConfig names the source regions and assigns each grid node to an
initial source via a source_mask vector (integer indices 1..n_sources).

ProvenanceState holds bed_mass_src[j, k, s] — the fraction of the surface-
layer mass at node j, fraction k that currently traces back to source s.

Invariant (checked in tests):
    sum_s( bed_mass_src[j, k, s] ) ≈ bed_mass[j, 1, k]   for all j, k

## Update algorithm

After each update_bed_composition! call the driver snapshots the before/after
bed_mass[:,1,:] change and attributes it via `step_provenance!`:

  * Net erosion at j (bed_mass decreased): remove proportionally from the
    existing source ledger at that node — the composition of eroded material
    mirrors whatever is in the active layer.
  * Net deposition at j (bed_mass increased): attribute new mass proportionally
    to the source composition at the "donor" node — whichever neighbor has the
    larger magnitude of net erosion (i.e. the upwind transport neighbor).
    Fallback: equal-weight average of both neighbors' source composition.

Only the surface layer (layer 1) is tracked.  This is sufficient for the
Hovmöller-style provenance visualisation and is the minimal defensible MVP.
==============================================================================#

"""
    ProvenanceConfig

Names the source regions and assigns the initial source index to every node.

Fields:
- `n_sources`     — number of distinct source regions
- `source_labels` — human-readable labels for each source
- `source_mask`   — integer vector of length n_x; value = source index (1..n_sources).
                    Must be allocated against the FULL grid (length config.grid.nn),
                    with entries beyond jmax left at any valid index (they are never read).
"""
struct ProvenanceConfig
    n_sources::Int
    source_labels::Vector{String}   # length n_sources
    source_mask::Vector{Int}        # length n_x; values in 1..n_sources
end

function ProvenanceConfig(source_labels::Vector{String}, source_mask::Vector{Int})
    n = length(source_labels)
    n >= 1 || throw(ArgumentError("n_sources must be ≥ 1"))
    all(1 .<= source_mask .<= n) ||
        throw(ArgumentError("source_mask values must be in 1..$(n)"))
    return ProvenanceConfig(n, source_labels, source_mask)
end

"""
    ProvenanceState

Mutable per-timestep provenance ledger.

- `bed_mass_src[j, k, s]` — kg/m², mass of surface-layer fraction k at node j
   that traces back to source s.  Layer index is omitted: only layer 1 (the
   active surface layer) is tracked.

Invariant: sum_s bed_mass_src[j,k,s] ≈ bed_mass[j,1,k]  ∀ j,k
"""
mutable struct ProvenanceState
    bed_mass_src::Array{Float64, 3}   # (n_x, n_frac, n_sources)
end

"""
    init_provenance(prov_config, state, config) -> ProvenanceState

Allocate and initialize ProvenanceState from the initial bed_mass.

At t=0, every node j receives 100% of its surface-layer mass from
source prov_config.source_mask[j].  All other source slots are zero.
"""
function init_provenance(prov_cfg::ProvenanceConfig,
                         state::CshoreState,
                         config::CshoreConfig)
    nn = config.grid.nn
    nf = nfractions(config.multifraction)
    ns = prov_cfg.n_sources
    jmax = state.jmax[1]   # use line-1 jmax as the domain extent

    bed_mass_src = zeros(Float64, nn, nf, ns)

    @inbounds for j in 1:jmax
        s = prov_cfg.source_mask[j]
        for k in 1:nf
            bed_mass_src[j, k, s] = state.bed_mass[j, 1, k]
        end
    end

    return ProvenanceState(bed_mass_src)
end

"""
    step_provenance!(prov, prov_cfg, bed_mass_before, state, config, l)

Update the provenance ledger after one morphodynamic sub-step.

Arguments:
- `prov`              — ProvenanceState (modified in place)
- `prov_cfg`          — ProvenanceConfig (read-only)
- `bed_mass_before`   — snapshot of state.bed_mass[:,1,:] BEFORE the sub-step
                         (shape: jmax × nf, only layer 1)
- `state`             — current CshoreState (after update_bed_composition!)
- `config`            — CshoreConfig
- `l`                 — cross-shore line index (currently only l=1 supported)

Algorithm
---------
For each node j and fraction k:

  dm = bed_mass_after[j,1,k] - bed_mass_before[j,k]

  If dm < 0 (erosion): remove from the source ledger proportionally to the
  current source composition at that node.

  If dm > 0 (deposition): attribute to the source composition of the
  "donor" neighbour.  The donor is whichever neighbour (j-1 or j+1) had
  the largest net erosion in this sub-step — a proxy for the upwind source.
  If both neighbours have equal (or zero) change, use an average of both.

After attribution, clamp all entries to ≥ 0 (numerical guard).
"""
function step_provenance!(prov::ProvenanceState,
                           prov_cfg::ProvenanceConfig,
                           bed_mass_before::Matrix{Float64},
                           state::CshoreState,
                           config::CshoreConfig,
                           l::Int)
    jmax_l = state.jmax[l]
    nf = size(state.bed_mass, 3)
    ns = prov_cfg.n_sources

    bms = prov.bed_mass_src   # (nn, nf, ns)

    # Net change in surface layer mass per node per fraction.
    # dm[j, k] = bed_mass_after[j,1,k] - bed_mass_before[j,k]
    # Use only layer 1.
    @inbounds for j in 2:(jmax_l - 1)
        for k in 1:nf
            dm = state.bed_mass[j, 1, k] - bed_mass_before[j, k]

            if dm < -1e-14
                # ---- Erosion: remove proportionally from source ledger ----
                total_src = 0.0
                for s in 1:ns
                    total_src += bms[j, k, s]
                end
                total_src_safe = max(total_src, 1e-14)
                for s in 1:ns
                    bms[j, k, s] += dm * (bms[j, k, s] / total_src_safe)
                end

            elseif dm > 1e-14
                # ---- Deposition: attribute to donor neighbour's composition ----
                # Pick the donor: whichever neighbour lost more mass (more erosion).
                # dm_left / dm_right = net change at j-1, j+1 (negative = erosion).
                dm_left  = state.bed_mass[j - 1, 1, k] - bed_mass_before[j - 1, k]
                dm_right = state.bed_mass[j + 1, 1, k] - bed_mass_before[j + 1, k]

                # Erosion magnitudes at the two neighbours
                ero_left  = dm_left  < 0 ? -dm_left  : 0.0
                ero_right = dm_right < 0 ? -dm_right : 0.0

                # Weighted mix: donor_frac[s] = fraction of deposited mass from source s
                w_left  = ero_left  / max(ero_left + ero_right, 1e-30)
                w_right = ero_right / max(ero_left + ero_right, 1e-30)

                # If neither neighbour is eroding, fall back to equal weight
                if ero_left + ero_right < 1e-30
                    w_left  = 0.5
                    w_right = 0.5
                end

                # Donor source fractions (normalised from their source ledger)
                jleft  = j - 1
                jright = j + 1

                total_left  = 0.0; total_right = 0.0
                for s in 1:ns
                    total_left  += bms[jleft,  k, s]
                    total_right += bms[jright, k, s]
                end
                total_left  = max(total_left,  1e-14)
                total_right = max(total_right, 1e-14)

                for s in 1:ns
                    frac_left  = bms[jleft,  k, s] / total_left
                    frac_right = bms[jright, k, s] / total_right
                    donor_frac = w_left * frac_left + w_right * frac_right
                    bms[j, k, s] += dm * donor_frac
                end
            end

            # Clamp negatives (numerical guard only — should be near-zero)
            for s in 1:ns
                if bms[j, k, s] < 0.0
                    bms[j, k, s] = 0.0
                end
            end
        end
    end
    return nothing
end

"""
    provenance_fractions(prov, j, k) -> Vector{Float64}

Return the normalized source fractions for surface-layer fraction k at node j.
Sum equals 1 (or all-zero if the node is empty).
"""
function provenance_fractions(prov::ProvenanceState, j::Int, k::Int)
    ns = size(prov.bed_mass_src, 3)
    total = 0.0
    for s in 1:ns
        total += prov.bed_mass_src[j, k, s]
    end
    fracs = Vector{Float64}(undef, ns)
    if total > 0
        for s in 1:ns
            fracs[s] = prov.bed_mass_src[j, k, s] / total
        end
    else
        fill!(fracs, 0.0)
    end
    return fracs
end
