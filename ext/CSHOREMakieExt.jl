#==============================================================================
CSHOREMakieExt — CairoMakie-backed implementations of CSHORE plotting stubs.

Loaded automatically when both `CSHORE` and `CairoMakie` are present in the
session (Julia 1.9+ package-extension mechanism). The corresponding stubs
live in `src/plotting.jl`; we add concrete methods here that take precedence
over the catch-all error fallback.

When CairoMakie is NOT loaded (e.g. inside the compiled `cshore-julia`
binary), this file is never read and the stubs raise a clear error pointing
the user at `using CairoMakie`.
==============================================================================#

module CSHOREMakieExt

import CSHORE
using CSHORE: CshoreState
import CairoMakie
import NCDatasets
import Dates
using Printf

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

const _PLOT_THEME = CairoMakie.Theme(
    backgroundcolor = :white,
    Axis = (
        xgridvisible  = true,
        ygridvisible  = true,
        xgridstyle    = :dash,
        ygridstyle    = :dash,
        xgridalpha    = 0.3,
        ygridalpha    = 0.3,
    ),
)

_with_theme(f) = CairoMakie.with_theme(f, _PLOT_THEME)

function _time_unit(total_seconds::Real)
    total_hours = total_seconds / 3600.0
    if total_hours < 200
        return 3600.0, "Time (hours)", "%.1f h"
    elseif total_hours < 2 * 365.25 * 24
        return 86400.0, "Time (days)", "%.1f d"
    else
        return 365.25 * 86400.0, "Time (years)", "%.2f yr"
    end
end

function _seconds_to_display(seconds_vec::AbstractVector{<:Real})
    total = length(seconds_vec) > 0 ? maximum(seconds_vec) - minimum(seconds_vec) : 0.0
    divisor, label, _ = _time_unit(total)
    return Float64.(seconds_vec) ./ divisor, label
end

function _format_time(seconds::Real, total_seconds::Real)
    divisor, _, _ = _time_unit(total_seconds)
    val = seconds / divisor
    if divisor == 3600.0
        return @sprintf("t = %.1f h", val)
    elseif divisor == 86400.0
        return @sprintf("t = %.1f d", val)
    else
        return @sprintf("t = %.2f yr", val)
    end
end

function _empty_figure(msg::AbstractString)
    @warn msg
    fig = CairoMakie.Figure(size = (400, 200))
    CairoMakie.Axis(fig[1,1]; title = msg)
    fig
end

function _safe_jmax(state::CshoreState, l::Int)
    l < 1 && return 0
    l > length(state.jmax) && return 0
    return state.jmax[l]
end

function _snapshot_indices(total::Int, n::Int)
    n = clamp(n, 1, total)
    n == 1 && return [total]
    n >= total && return collect(1:total)
    return unique(round.(Int, range(1, total; length = n)))
end

# ---------------------------------------------------------------------------
# 1. plot_profile
# ---------------------------------------------------------------------------

function CSHORE.plot_profile(x::AbstractVector, z::AbstractVector;
                             hardbottom_z::Union{Nothing,AbstractVector}=nothing,
                             title::AbstractString="",
                             swl::Real=0.0)
    length(x) == 0 && return _empty_figure("plot_profile: empty input")
    _with_theme() do
        fig = CairoMakie.Figure(size=(900, 400))
        ax  = CairoMakie.Axis(fig[1,1];
            xlabel = "Cross-shore distance (m)",
            ylabel = "Elevation (m)",
            title  = isempty(title) ? "Cross-shore profile" : title,
        )
        if hardbottom_z !== nothing && length(hardbottom_z) == length(x)
            CairoMakie.band!(ax, x, hardbottom_z, z;
                color = (:sandybrown, 0.4), label = "Erodible sediment")
            CairoMakie.lines!(ax, x, hardbottom_z;
                color = :sienna, linewidth = 1.5, linestyle = :dash,
                label = "Hardbottom")
        end
        CairoMakie.lines!(ax, x, z; color = :black, linewidth = 2.0, label = "Bed")
        CairoMakie.hlines!(ax, [swl]; color = :dodgerblue, linewidth = 1.0,
            linestyle = :dashdot, label = "SWL = $(round(swl; digits=2)) m")
        z_clip = [min(zi, swl) for zi in z]
        CairoMakie.band!(ax, x, z_clip, fill(swl, length(x));
            color = (:dodgerblue, 0.15))
        CairoMakie.axislegend(ax; position = :rt, framevisible = false)
        fig
    end
end

# ---------------------------------------------------------------------------
# 2. plot_profile_evolution
# ---------------------------------------------------------------------------

function CSHORE.plot_profile_evolution(nc_file::AbstractString;
                                       snapshots::Int=5,
                                       title::AbstractString="")
    isfile(nc_file) || error("NetCDF file not found: $nc_file")
    _with_theme() do
        NCDatasets.NCDataset(nc_file, "r") do ds
            x    = Array(ds["x"])
            zb_raw = Array(ds["zb"])
            zb   = permutedims(zb_raw)
            time_raw = Array(ds["time"])
            time = if eltype(time_raw) <: Dates.AbstractDateTime
                t0 = time_raw[1]
                Float64[Dates.value(t - t0) / 1000.0 for t in time_raw]
            else
                Float64.(time_raw)
            end
            nt   = size(zb, 1)
            nt < 1 && return _empty_figure("plot_profile_evolution: no time steps")

            idxs = _snapshot_indices(nt, snapshots)
            cmap = CairoMakie.cgrad(:viridis, length(idxs); categorical = true)

            fig = CairoMakie.Figure(size = (900, 450))
            ax  = CairoMakie.Axis(fig[1,1];
                xlabel = "Cross-shore distance (m)",
                ylabel = "Elevation (m)",
                title  = isempty(title) ? "Profile evolution" : title,
            )

            total_s = time[end] - time[1]
            for (ci, ti) in enumerate(idxs)
                lbl = _format_time(time[ti], total_s)
                CairoMakie.lines!(ax, x, view(zb, ti, :);
                    color = cmap[ci], linewidth = 1.5, label = lbl)
            end
            CairoMakie.axislegend(ax; position = :rt, framevisible = false)
            fig
        end
    end
end

# ---------------------------------------------------------------------------
# 3. plot_wave_field
# ---------------------------------------------------------------------------

function CSHORE.plot_wave_field(state::CshoreState; l::Int=1)
    jm = _safe_jmax(state, l)
    jm < 2 && return _empty_figure("plot_wave_field: jmax < 2")
    _with_theme() do
        x = view(state.xb, 1:jm)
        fig = CairoMakie.Figure(size = (900, 700))

        ax1 = CairoMakie.Axis(fig[1,1]; ylabel = "Hrms (m)", title = "Wave field (line $l)")
        CairoMakie.lines!(ax1, x, view(state.hrms, 1:jm); color = :teal, linewidth = 1.5)

        ax2 = CairoMakie.Axis(fig[2,1]; ylabel = "Setup (m)")
        CairoMakie.lines!(ax2, x, view(state.wsetup, 1:jm); color = :coral, linewidth = 1.5)

        ax3 = CairoMakie.Axis(fig[3,1]; xlabel = "Cross-shore distance (m)",
                                         ylabel = "Depth (m)")
        CairoMakie.lines!(ax3, x, view(state.h, 1:jm); color = :steelblue, linewidth = 1.5)

        CairoMakie.linkxaxes!(ax1, ax2, ax3)
        CairoMakie.hidexdecorations!(ax1; grid = false)
        CairoMakie.hidexdecorations!(ax2; grid = false)
        fig
    end
end

# ---------------------------------------------------------------------------
# 4. plot_transport
# ---------------------------------------------------------------------------

function CSHORE.plot_transport(state::CshoreState; l::Int=1, fraction::Int=1)
    jm = _safe_jmax(state, l)
    jm < 2 && return _empty_figure("plot_transport: jmax < 2")
    _with_theme() do
        x   = view(state.xb, 1:jm)
        fig = CairoMakie.Figure(size = (900, 550))

        ax1 = CairoMakie.Axis(fig[1,1]; ylabel = "Transport rate (m³/m/s)",
                                         title  = "Sediment transport (fraction $fraction, line $l)")
        CairoMakie.lines!(ax1, x, view(state.qbx, 1:jm, fraction);
            color = :darkorange, linewidth = 1.5, label = "QBX (bedload)")
        CairoMakie.lines!(ax1, x, view(state.qsx, 1:jm, fraction);
            color = :purple, linewidth = 1.5, label = "QSX (suspended)")
        CairoMakie.axislegend(ax1; position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(fig[2,1]; xlabel = "Cross-shore distance (m)",
                                         ylabel = "Transport rate (m³/m/s)")
        CairoMakie.lines!(ax2, x, view(state.qby, 1:jm, fraction);
            color = :darkorange, linewidth = 1.5, linestyle = :dash, label = "QBY (bedload)")
        CairoMakie.lines!(ax2, x, view(state.qsy, 1:jm, fraction);
            color = :purple, linewidth = 1.5, linestyle = :dash, label = "QSY (suspended)")
        CairoMakie.axislegend(ax2; position = :rt, framevisible = false)

        CairoMakie.linkxaxes!(ax1, ax2)
        CairoMakie.hidexdecorations!(ax1; grid = false)
        fig
    end
end

# ---------------------------------------------------------------------------
# 5. plot_hovmoller
# ---------------------------------------------------------------------------

function CSHORE.plot_hovmoller(nc_file::AbstractString;
                               var::Union{Symbol,String}=:zb,
                               title::AbstractString="")
    isfile(nc_file) || error("NetCDF file not found: $nc_file")
    varname = string(var)
    _with_theme() do
        NCDatasets.NCDataset(nc_file, "r") do ds
            haskey(ds, varname) || error("Variable '$varname' not found in $nc_file")
            x    = Array(ds["x"])
            time_raw = Array(ds["time"])
            time_s = if eltype(time_raw) <: Dates.AbstractDateTime
                t0 = time_raw[1]
                Float64[Dates.value(t - t0) / 1000.0 for t in time_raw]
            else
                Float64.(time_raw)
            end
            time_display, time_label = _seconds_to_display(time_s)
            data_raw = Array(ds[varname])
            data = permutedims(data_raw)

            is_d50 = occursin("d50", varname)
            display_data = is_d50 ? data .* 1e3 : data
            cb_label = is_d50 ? "$varname (mm)" : varname

            size(data, 1) < 1 && return _empty_figure("plot_hovmoller: no time steps")

            fig = CairoMakie.Figure(size = (900, 500))
            ax  = CairoMakie.Axis(fig[1,1];
                xlabel = "Cross-shore distance (m)",
                ylabel = time_label,
                title  = isempty(title) ? varname : title,
            )
            hm = CairoMakie.heatmap!(ax, x, time_display, permutedims(display_data);
                colormap = is_d50 ? CairoMakie.Reverse(:RdYlBu) : :balance)
            CairoMakie.Colorbar(fig[1,2], hm; label = cb_label)
            fig
        end
    end
end

# ---------------------------------------------------------------------------
# 6. plot_mass_balance
# ---------------------------------------------------------------------------

function CSHORE.plot_mass_balance(nc_file::AbstractString)
    isfile(nc_file) || error("NetCDF file not found: $nc_file")
    _with_theme() do
        NCDatasets.NCDataset(nc_file, "r") do ds
            haskey(ds, "bed_mass") || error("'bed_mass' not found in $nc_file")
            bm   = Array(ds["bed_mass"])
            time_raw = Array(ds["time"])
            time_s = if eltype(time_raw) <: Dates.AbstractDateTime
                t0 = time_raw[1]
                Float64[Dates.value(t - t0) / 1000.0 for t in time_raw]
            else
                Float64.(time_raw)
            end
            time_display, time_label = _seconds_to_display(time_s)

            grain_labels = if haskey(ds, "fraction")
                gs = Array(ds["fraction"])
                [string(round(g * 1e3, digits=3), " mm") for g in gs]
            else
                nothing
            end

            nx_bm, nlay, nf, nt = size(bm)
            nt < 1 && return _empty_figure("plot_mass_balance: no time steps")

            total = zeros(nt, nf)
            for ti in 1:nt, k in 1:nf
                for lay in 1:nlay, j in 1:nx_bm
                    total[ti, k] += max(0.0, bm[j, lay, k, ti])
                end
            end

            fig = CairoMakie.Figure(size = (900, 400))
            ax  = CairoMakie.Axis(fig[1,1];
                xlabel = time_label,
                ylabel = "Total bed mass (kg/m)",
                title  = "Per-fraction mass balance",
            )
            for k in 1:nf
                lbl = grain_labels !== nothing ? grain_labels[k] : "Fraction $k"
                CairoMakie.lines!(ax, time_display, view(total, :, k);
                    linewidth = 1.5, label = lbl)
            end
            CairoMakie.axislegend(ax; position = :rt, framevisible = false)
            fig
        end
    end
end

# ---------------------------------------------------------------------------
# 7. plot_swash
# ---------------------------------------------------------------------------

function CSHORE.plot_swash(state::CshoreState; l::Int=1)
    jm = _safe_jmax(state, l)
    j1 = max(1, state.jwd)
    j2 = min(jm, state.jr)
    j2 <= j1 && return _empty_figure("plot_swash: no swash zone (jwd=$j1, jr=$j2)")
    _with_theme() do
        rng = j1:j2
        x   = view(state.xb, rng)
        fig = CairoMakie.Figure(size = (900, 500))

        ax1 = CairoMakie.Axis(fig[1,1]; ylabel = "Wet probability",
                                         title  = "Swash zone diagnostics (line $l)")
        CairoMakie.lines!(ax1, x, view(state.pwet, rng);
            color = :teal, linewidth = 1.5)
        CairoMakie.ylims!(ax1, -0.05, 1.05)

        ax2 = CairoMakie.Axis(fig[2,1]; xlabel = "Cross-shore distance (m)",
                                         ylabel = "Swash depth (m)")
        CairoMakie.lines!(ax2, x, view(state.hwd, rng);
            color = :steelblue, linewidth = 1.5)

        CairoMakie.linkxaxes!(ax1, ax2)
        CairoMakie.hidexdecorations!(ax1; grid = false)
        fig
    end
end

# ---------------------------------------------------------------------------
# 8. plot_thermal
# ---------------------------------------------------------------------------

function CSHORE.plot_thermal(nc_file::AbstractString; title::AbstractString="")
    isfile(nc_file) || error("NetCDF file not found: $nc_file")
    _with_theme() do
        NCDatasets.NCDataset(nc_file, "r") do ds
            has_alt = haskey(ds, "ALT")
            has_alt || return _empty_figure("plot_thermal: no ALT variable — run with THERMAL=true")

            x    = Array(ds["x"])
            time_raw = Array(ds["time"])
            time_s = if eltype(time_raw) <: Dates.AbstractDateTime
                t0 = time_raw[1]
                Float64[Dates.value(t - t0) / 1000.0 for t in time_raw]
            else
                Float64.(time_raw)
            end
            time_disp, time_label = _seconds_to_display(time_s)

            alt_nc  = permutedims(Array(ds["ALT"]))
            zbh_nc  = permutedims(Array(ds["zb_hard"]))
            zb_nc   = permutedims(Array(ds["zb"]))

            has_tsurf = haskey(ds, "T_surface")
            tsurf_nc = has_tsurf ? permutedims(Array(ds["T_surface"])) : nothing

            nt, nx = size(alt_nc)
            ttl = isempty(title) ? "Thermal / Permafrost Diagnostics" : title

            fig = CairoMakie.Figure(size=(1100, 900))

            ax1 = CairoMakie.Axis(fig[1,1]; xlabel="Cross-shore distance (m)",
                ylabel=time_label, title="Active Layer Thickness (m)")
            hm1 = CairoMakie.heatmap!(ax1, x, time_disp, permutedims(alt_nc);
                colormap=:YlOrRd)
            CairoMakie.Colorbar(fig[1,2], hm1; label="ALT (m)")

            if has_tsurf
                ax2 = CairoMakie.Axis(fig[2,1]; xlabel="Cross-shore distance (m)",
                    ylabel=time_label, title="Surface Temperature (°C)")
                hm2 = CairoMakie.heatmap!(ax2, x, time_disp, permutedims(tsurf_nc);
                    colormap=:RdBu)
                CairoMakie.Colorbar(fig[2,2], hm2; label="T (°C)")
            end

            ax3 = CairoMakie.Axis(fig[3,1]; xlabel="Cross-shore distance (m)",
                ylabel="Elevation (m)", title="Final Profile + Permafrost Table")
            CairoMakie.band!(ax3, x, fill(minimum(zbh_nc[end,:]) - 1, nx), zbh_nc[end, :],
                color=(:lightblue, 0.4), label="Frozen ground")
            CairoMakie.band!(ax3, x, zbh_nc[end, :], zb_nc[end, :],
                color=(:tan, 0.5), label="Active layer")
            CairoMakie.lines!(ax3, x, zb_nc[end, :], color=:black, linewidth=1.5, label="Bed surface")
            CairoMakie.lines!(ax3, x, zbh_nc[end, :], color=:steelblue, linewidth=1.2,
                linestyle=:dash, label="Permafrost table")
            CairoMakie.axislegend(ax3, position=:lt)

            ax4 = CairoMakie.Axis(fig[4,1]; xlabel=time_label,
                ylabel="Max ALT (m)", title="Maximum Active Layer Thickness")
            alt_max = [maximum(alt_nc[ti, :]) for ti in 1:nt]
            CairoMakie.lines!(ax4, time_disp, alt_max, color=:firebrick, linewidth=1.5)
            CairoMakie.band!(ax4, time_disp, fill(0.0, nt), alt_max,
                color=(:orange, 0.2))

            CairoMakie.Label(fig[0, :], ttl, fontsize=15, font=:bold)
            fig
        end
    end
end

# ---------------------------------------------------------------------------
# 9. save_figure
# ---------------------------------------------------------------------------

function CSHORE.save_figure(fig::CairoMakie.Figure, path::AbstractString; kwargs...)
    CairoMakie.save(path, fig; kwargs...)
    @info "Figure saved to $path"
    path
end

end # module CSHOREMakieExt
