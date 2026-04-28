using CairoMakie
using LaTeXStrings
using CSV, DataFrames
using HDF5
using Printf
using Glob
using Statistics


loadtxt(str) = Matrix(CSV.read(str, DataFrame; delim=' ', header=false, ignorerepeated=true))
find_nearest(array :: AbstractVector, value :: Float64) = argmin(abs.(array .- value))
extract_val(prefix, str) = parse(Float64, match(Regex("$(prefix)_(-?\\d*\\.?\\d+)"), str).captures[1])

pic_dir = joinpath(pwd(), "pic_dir")
mkpath(pic_dir)
my_colors = cgrad(:glasbey_bw_minc_20_maxl_70_n256, 256, categorical=true)

custom_theme = Theme(
    fontsize = 24,
    Axis = (
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xlabelsize = 24,
        ylabelsize = 24,
        xtickformat = x -> [@sprintf("%.1f", val) for val in x],
        markersize = 10,
        linewidth = 1.75,
    ),
    palette = (
        color = my_colors,
    ),
    Legend = (
        titlesize = 20,
        labelsize = 20,
        markersize = 10,
        framevisible = false,
        tellwidth = false,
        tellheight = false,
    ),
)

set_theme!(merge(theme_latexfonts(), custom_theme))


## common functions

@recipe(CurvedArrow) do scene
    Theme(
        color = :black,
    )
end

function Makie.plot!(p::CurvedArrow)
    scene = Makie.get_scene(p)

    points = lift(p, scene.camera.projectionview, p.model, Makie.transform_func(p),
          scene.viewport, p[1], p[2]) do _, _, _, _, p1, p2

        return Makie.project.(Ref(scene), (p1, p2))
    end

    arc = lift(points) do (p1, p2)
        len = Makie.norm(p2 - p1)
        EllipticalArc(
            p1...,
            p2...,
            len/1.3,
            len/1.3,
            0,
            false,
            false,
        )
    end

    path = lift(arc) do arc
        BezierPath([
            MoveTo(points[][1]),
            arc,
        ])
    end

    trimarker = BezierPath([MoveTo(0, 0), LineTo(0.5, -1), LineTo(-0.5, -1), ClosePath()])

    markerangle = lift(arc) do arc
        arc.a2 + pi
    end

    markerangle_start = lift(arc) do arc
        arc.a1
    end

    lines!(p, path, space = :pixel, color = p.color, linestyle = Linestyle([0.0, 4.0, 8.0, 9.5]))
    scatter!(p, p[2], marker = trimarker, rotation = markerangle, color = p.color)

    scatter!(p, p[1],
    marker = trimarker,
    rotation = markerangle_start,
    color = p.color,
    )
end

function plot_electrode_diagram!(ax; levels=(-1.0, 1.0), voltage=2.0, lw=6,
    electrodeY=(-3, 3), connect_both=false, if_oscillator=false,
    draw_battery=true, battery_drop=1.5,
    battery_gap=0.02, battery_short=0.5, battery_long=1.5,
    battery_wire_lw=1.5, battery_lw=1.5)

    X = (-0.2, 0.2)
    Y = electrodeY
    aspect = 0.25
    mus = (voltage, -voltage)
    bar_width = 0.15
    electrode_fill_color = :navyblue
    electrode_empty_color = (:navyblue, 0.2)

    empty!(ax)
    bez(p0, p1, p2, t) = @. (1 - t)^2 * p0 + 2 * (1 - t) * t * p1 + t^2 * p2
    ang(p1, p2) = atan((p2[2] - p1[2]) * aspect, p2[1] - p1[1])
    gaussian(y, amp, sigma) = amp * exp(-0.5 * (y / sigma)^2)

    gap_x = range(X[1], X[2], length=100)
    sigma_h = 0.08

    x_offsets = [0.05, -0.05]
    for (i, (lvl, color)) in enumerate(zip(levels, [:red, :darkgreen]))
        x_shifted = ([X...] .* 0.4) .+ x_offsets[i]
        lines!(ax, x_shifted, [lvl, lvl]; color=color, linewidth=lw)
    end

    mid_lvl = sum(levels) / 2
    gap_x = range(X[1], X[2], length=100)
    sigma_h = 0.073
    amp = 3.0
    y_top = [amp .* exp(-0.5 * (x / sigma_h)^2) + mid_lvl for x in gap_x]
    y_bot = [-amp .* exp(-0.5 * (x / sigma_h)^2) + mid_lvl for x in gap_x]

    fill_points = [
        [Point2f(x, yt) for (x, yt) in zip(gap_x, y_top)]...,
        [Point2f(x, yb) for (x, yb) in zip(reverse(gap_x), reverse(y_bot))]...
    ]

    poly!(ax, fill_points, color=(:darkgoldenrod1, 0.1), strokewidth=0)
    lines!(ax, gap_x, y_top, color=(:darkgoldenrod1, 0.3), linewidth=1.5)
    lines!(ax, gap_x, y_bot, color=(:darkgoldenrod1, 0.3), linewidth=1.5)

    for (i, (lvl, color)) in enumerate(zip(levels, [:red, :darkgreen]))
        x_shifted = ([X...] .* 0.4) .+ x_offsets[i]
        lines!(ax, x_shifted, [lvl, lvl]; color=color, linewidth=lw)

        if if_oscillator
            line_mid = mean(x_shifted)
            if i == 1
                x_vals = range(line_mid - 0.08, line_mid + 0.08, length=50)
                y_vals = 350 .* (x_vals .- line_mid).^2 .+ lvl
                lines!(ax, x_vals, y_vals, linewidth=2, color="brown")
            elseif i == 2
                x_vals = range(line_mid - 0.015, line_mid + 0.1, length=60)
                De, a = 2.0, 50.0
                y_vals = De .* (1 .- exp.(-a .* (x_vals .- line_mid))).^2 .+ lvl
                lines!(ax, x_vals, y_vals, linewidth=2, color="brown")
            end
        end
    end

    for (i, x) in enumerate(X)
        mu = mus[i]
        sgn = (i == 1) ? -1 : 1
        x_outer = x + sgn * bar_width

        poly!(ax, [Point2f(x, Y[1]), Point2f(x_outer, Y[1]),
                Point2f(x_outer, mu), Point2f(x, mu)],
            color=electrode_fill_color, strokecolor=electrode_fill_color,
            strokewidth=1.5)

        poly!(ax, [Point2f(x, mu), Point2f(x_outer, mu),
                Point2f(x_outer, Y[2]), Point2f(x, Y[2])],
            color=electrode_empty_color, strokecolor=electrode_empty_color,
            strokewidth=0.0)

        mu_label = (i == 1) ? L"\mu_R" : L"\mu_L"
        text!(ax, Point2f(x + sgn * bar_width / 2, minimum(Y) + 0.25),
            text = mu_label,
            align = (:center, :below),
            fontsize = 25,
            color = :white)

        target_indices = connect_both ? [1, 2] : [(3 - i)]

        for current_target_idx in target_indices
            lvl = levels[current_target_idx]

            line_mid = mean(([X...] .* 0.4) .+ x_offsets[current_target_idx])

            p_start = Point2f(x, mu)
            p_end   = Point2f(line_mid, lvl)

            if current_target_idx == 1
                curvedarrow!(ax, p_start, p_end; color=:black)
            else
                curvedarrow!(ax, p_end, p_start; color=:black)
            end
        end
    end

    if draw_battery
        y_wire = Y[1] - battery_drop

        left_bottom_mid_x  = X[1] - bar_width / 2
        right_bottom_mid_x = X[2] + bar_width / 2

        x_bar1 = - battery_gap / 2
        x_bar2 = + battery_gap / 2

        x_left_stop   = x_bar1
        x_right_start = x_bar2

        lines!(ax,
            [left_bottom_mid_x, left_bottom_mid_x, x_left_stop],
            [Y[1], y_wire, y_wire];
            color = :black, linewidth = battery_wire_lw
        )

        lines!(ax,
            [x_bar1, x_bar1],
            [y_wire - battery_long/2, y_wire + battery_long/2];
            color = :black, linewidth = battery_lw
        )
        lines!(ax,
            [x_bar2, x_bar2],
            [y_wire - battery_short/2, y_wire + battery_short/2];
            color = :black, linewidth = battery_lw
        )

        lines!(ax,
            [x_right_start, right_bottom_mid_x, right_bottom_mid_x],
            [y_wire, y_wire, Y[1]];
            color = :black, linewidth = battery_wire_lw
        )
    end

    hidedecorations!(ax)
    hidespines!(ax)
    ylims!(ax, Y[1] - (draw_battery ? battery_drop + battery_long/2 + 0.1 : 0), Y[2])
    xlims!(ax, X[1] - bar_width, X[2] + bar_width)
    tightlimits!(ax)

    return ax
end

function get_steady_state(dir_path, is_heom=false)
    file = is_heom ? "current.HEOM.txt" : "current_average.txt"
    filepath = joinpath(dir_path, file)
    isfile(filepath) || return NaN

    f1 = loadtxt(filepath)
    times = f1[:, 1]
    current = -0.5 .* (f1[:, 2] .- f1[:, 3])
    start_idx = find_nearest(times, times[end] - 3.0)
    return mean(current[start_idx:end])
end


function plot_steady_state!(ax1; set="SetII", target_dir="scan_mu_nuclear",
    plot_x="lamda", plot_every=1, left_mu_limit = -5.0, colors = my_colors,
    if_colorbar=false, linestyle=nothing)

    target_dir = isabspath(target_dir) ? target_dir : joinpath(pwd(), "AA", target_dir)
    set_dir = joinpath(target_dir, set)
    lamda_dirs = sort(glob("lamda*", set_dir), by = x -> extract_val("lamda", x))

    lamda_vals = [extract_val("lamda", l) for l in lamda_dirs]
    mu_vals = zeros(length(glob("mu*", lamda_dirs[1])))
    steady_state_curr = zeros((length(lamda_dirs), length(mu_vals)))

    for (l, lamda) in enumerate(lamda_dirs)
        mu_dirs = sort(glob("mu*", lamda), by = x -> extract_val("mu", x))
        for (m, mu) in enumerate(mu_dirs)
            is_heom = contains(lowercase(mu), "heom")
            mu_val = extract_val("mu", mu)
            mu_vals[m] = is_heom ? -mu_val : mu_val
            steady_state_curr[l, m] = get_steady_state(mu, is_heom)
        end
    end

    k = 0
    lamda_idx = 1:plot_every:length(lamda_vals)

    if plot_x == "mu"
        for i in eachindex(mu_vals)
            if 0.0 < mu_vals[i] <= 7.0
                k += 1
                valid_mu_len = length(findall(x -> (0.0 < x <= 7.0), mu_vals))
                label = latexstring("\\mu_R/\\Gamma = $( @sprintf("%.3f", mu_vals[i]) )")
                color_idx = isnothing(linestyle) ? colors[k] : colors[valid_mu_len + 1 - k]

                if isnothing(linestyle)
                    scatterlines!(ax1, lamda_vals[lamda_idx], steady_state_curr[lamda_idx, i],
                         label=label, color=color_idx, markersize=12)
                else
                    lines!(ax1, lamda_vals[lamda_idx], steady_state_curr[lamda_idx, i],
                         label=label, color=color_idx, linestyle=linestyle)
                end
            end
        end
    elseif plot_x == "lamda"
        valid_mu_idx = findall(x -> (x >= left_mu_limit && x <= 7.0), mu_vals)
        for i in lamda_idx
            label = latexstring("\\lambda/\\Gamma = $( @sprintf("%.3f", lamda_vals[i]) )")

            if isnothing(linestyle)
                scatterlines!(ax1, mu_vals[valid_mu_idx], steady_state_curr[i, valid_mu_idx],
                     label=label, color=colors[i], markersize=12)
            else
                lines!(ax1, mu_vals[valid_mu_idx], steady_state_curr[i, valid_mu_idx],
                     label=label, color=colors[i], linestyle=linestyle)
            end
        end
    end

    return ax1, plot_x == "lamda" ? lamda_vals : mu_vals[findall(x -> (0.0 < x <= 7.0), mu_vals)]
end


function plot_color_bar(ax, lamda_vals;
        loc=(0.8, 0.15), myticks=nothing, label = L"\mu_R/\Gamma",
        tickformat=nothing, colors = colors, width = 0.5)
    kwargs = (
        bbox = ax.scene.viewport,
        colormap = colors[1:length(lamda_vals)],
        limits = isnothing(myticks) ? extrema(lamda_vals) : extrema(myticks),
        ticklabelsize = 18,
        vertical = false,
        halign = loc[1],
        valign = loc[2],
        width = Relative(width),
        label = label,
        labelpadding = -70,
    )

    if isnothing(myticks)
        if isnothing(tickformat)
            return Colorbar(fig; kwargs...)
        else
            return Colorbar(fig; kwargs..., tickformat=tickformat)
        end
    else
        if isnothing(tickformat)
            return Colorbar(fig; kwargs..., ticks=myticks)
        else
            return Colorbar(fig; kwargs..., tickformat=tickformat, ticks=myticks)
        end
    end
end

## figure 1
begin
    delta = 1.0
    omega_t = 0.1 * delta; omega_c = 0.2 * delta

    epsilon1 = -delta/2; epsilon2 = delta/2

    kappa1 = -0.2 * delta; kappa2 = 0.3 * delta
    gamma = 0.05 * delta

    function adiabatic_potentials(Qt, Qc)
        V0 = 0.5 * omega_t * abs2(Qt) + 0.5 * omega_c * abs2(Qc)

        V11 = V0 + epsilon1 + kappa1 * Qt
        V22 = V0 + epsilon2 + kappa2 * Qt
        V12 = gamma * Qc

        mean_V = (V11 + V22) * 0.5
        diff_V = (V11 - V22) * 0.5
        gap = sqrt(abs2(diff_V) + abs2(V12))
        return mean_V - gap, mean_V + gap
    end

    grid_size = 101
    Qt_range = range(-4.0, 0.0, length=grid_size)
    Qc_range = range(-2.0, 2.0, length=grid_size)

    E_minus = Array{Float64}(undef, (grid_size, grid_size))
    E_plus  = zeros(Float64, grid_size, grid_size)

    for (i, qt) in enumerate(Qt_range)
        for (j, qc) in enumerate(Qc_range)
            pot = adiabatic_potentials(qt, qc)
            E_minus[i, j] = pot[1]
            E_plus[i, j]  = pot[2]
        end
    end

    top_col_width = 0.8

    fig = Figure(size=(1500, 450), figure_padding=(10, 10, 10, 0))

    ax1 = Axis3(fig[1:2, 1],
            xlabel = L"Q_t",
            ylabel = L"Q_c",
            zlabel = "",
            azimuth = deg2rad(45),
            elevation = deg2rad(15),
            aspect=(2,2,2),
            xticks=LinearTicks(3),
            yticks=LinearTicks(3),
            title=L"\textrm{Adiabatic Potential,}\, \mathbf{h}(Q_c,Q_t)",
            titlegap=-25,
            )

    surface!(ax1, Qt_range, Qc_range, E_minus, colormap=:davos, alpha=0.8, rasterize = 5)
    surface!(ax1, Qt_range, Qc_range, E_plus, colormap=:acton, alpha=0.8, rasterize = 5)

    g12 = GridLayout(fig[1, 2])
    colsize!(g12, 1, Relative(top_col_width))
    ax2 = Axis(g12[1, 1])
    # ax2 = Axis(fig[1,2])
    plot_electrode_diagram!(ax2; voltage=1.8)

    ax3 = Axis(fig[2,2], ylabel=L"I_{SS}", xlabel=L"\mu_R/\Gamma", ylabelrotation=0)
    colors = cgrad(:coolwarm, 9, categorical = true)
    _, lamda_vals = plot_steady_state!(ax3; plot_x ="lamda",
        target_dir="scan_mu",
        set="SetIII",
        colors=colors,
        left_mu_limit=-3.0)
    plot_color_bar(ax3, lamda_vals; loc=(0.1, 0.85))

    g13 = GridLayout(fig[1, 3])
    colsize!(g13, 1, Relative(top_col_width))
    ax4 = Axis(g13[1, 1])
    # ax4 = Axis(fig[1,3])
    plot_electrode_diagram!(ax4; voltage=1.4, if_oscillator=true)

    ax5 = Axis(fig[2,3], xlabel=L"\mu_R/\Gamma")
    colors = cgrad(:coolwarm, 9, categorical = true)
    _, lamda_vals = plot_steady_state!(ax5; plot_x ="lamda",
        set="SetIII", colors=colors,
        left_mu_limit=-3.0,
        target_dir="scan_mu_nuclear",)
    # plot_color_bar(ax5, lamda_vals; loc=(0.1, 0.85))


    # ax6 = Axis(fig[1,4])
    g14 = GridLayout(fig[1, 4])
    colsize!(g14, 1, Relative(top_col_width))
    ax6 = Axis(g14[1, 1])
    plot_electrode_diagram!(ax6; voltage=0.8, if_oscillator=true, connect_both=true)

    ax7 = Axis(fig[2,4], xlabel=L"\mu_R/\Gamma")
    colors = cgrad(:coolwarm, 9, categorical = true)
    _, lamda_vals = plot_steady_state!(ax7; plot_x ="lamda",
        target_dir="scan_mu_nuclear_check_coupled",
        set="SetIII",
        plot_every=4, colors=colors,
        )
    # plot_color_bar(ax7, lamda_vals; loc=(0.1, 0.85))
    text!(ax7, 0.6, 0.2; text="Set III", space=:relative, font = :bold,)

    Label(fig[1, 1, TopLeft()], "(A)", font=:bold, padding=(0, 5, -50, 0), halign=:right)
    Label(fig[1, 2, TopLeft()], "(B)", font=:bold, padding=(0, 5, -50, 0), halign=:right)
    Label(fig[1, 3, TopLeft()], "(C)", font=:bold, padding=(0, 5, -50, 0), halign=:right)
    Label(fig[1, 4, TopLeft()], "(D)", font=:bold, padding=(0, 5, -50, 0), halign=:right)

    rowsize!(fig.layout, 1, Auto(0.5))
    rowsize!(fig.layout, 2, Auto(1.0))
    rowgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 1, 0)

    save(joinpath(pic_dir, "nuclear_pot_system_compare.pdf"), fig)
    fig
end


## figure 2

function plot_nuclear_modes_transformed(; sampling="AA", target_dir="scan_mu_nuclear",
        set="SetII", lamda_dir="lamda_0.0750")

    data_dir = joinpath(pwd(), sampling, target_dir, set, lamda_dir)
    mu_dirs = glob("mu_*", data_dir)
    mu_vals = [parse(Float64, match(r"mu_(-?\d+\.?\d*)", p).captures[1]) for p in mu_dirs]
    filted_mu_vals = [x for x in mu_vals if 0.8 <= x <= 3.0]
    filted_mu_dirs = [joinpath(data_dir, @sprintf("mu_%.1f", x)) for x in filted_mu_vals]

    raw_qc_data = Dict{Int,Vector{Float64}}()
    raw_qt_data = Dict{Int,Vector{Float64}}()

    for (mu_idx, mu_vals) in enumerate(filted_mu_vals)
        datafile = joinpath(filted_mu_dirs[mu_idx], "qc_qt_average.h5")
        @assert isfile(datafile)
        h5open(datafile, "r") do f
            raw = read(f["data"])
            raw_qc_data[mu_idx] = raw[:, 1]
            raw_qt_data[mu_idx] = raw[:, 2]
        end
    end

    fig = Figure(size=(750,350))

    ax1 = Axis3(fig[1,1], xlabel=L"Q_c", ylabel=L"\mu_R/\Gamma", zlabel="",
                yticks=LinearTicks(5),
                azimuth = deg2rad(-49),
                elevation = deg2rad(25),
                aspect=(3, 3, 1.5),
                zticks=LinearTicks(4))
    density_color_qc = :orchid4
    # density_color_qc = :black
    for (i, mu) in enumerate(filted_mu_vals)
        plt = density!(ax1, raw_qc_data[i], strokewidth=1.2,
                       color=(density_color_qc, 0.2),
                       strokecolor=density_color_qc)
        rotate!(plt, Vec3f(1, 0, 0), pi/2)
        translate!(plt, 0, mu, 0)
    end

    qc_means = [mean(raw_qc_data[i]) for i in eachindex(filted_mu_vals)]
    qc_stds  = [std(raw_qc_data[i]) for i in eachindex(filted_mu_vals)]

    ylims!(ax1, 0.5, 3.2)
    zlims!(ax1, 0.0, 0.6)

    ax2 = Axis3(fig[1,2], xlabel=L"Q_t", ylabel=L"\mu_R/\Gamma", zlabel="",
                yticks=LinearTicks(5),
                azimuth = deg2rad(-49),
                elevation = deg2rad(25),
                aspect=(3, 3, 1.5),
                )
    density_color_qt = :turquoise4
    # density_color_qt = :black
    for (i, mu) in enumerate(filted_mu_vals)
        plt = density!(ax2, raw_qt_data[i], strokewidth=1.2,
                       color=(density_color_qt, 0.2),
                       strokecolor=density_color_qt)
        rotate!(plt, Vec3f(1, 0, 0), pi/2)
        translate!(plt, 0, mu, 0)
    end

    qt_means = [mean(raw_qt_data[i]) for i in eachindex(filted_mu_vals)]
    qt_stds  = [std(raw_qt_data[i]) for i in eachindex(filted_mu_vals)]

    z_floor = 0.001
    lower_pts = [Point3f(m - s, y, z_floor) for (m, s, y) in zip(qt_means, qt_stds, filted_mu_vals)]
    upper_pts = [Point3f(m + s, y, z_floor) for (m, s, y) in zip(qt_means, qt_stds, filted_mu_vals)]
    mean_pts  = [Point3f(m, y, z_floor) for (m, y) in zip(qt_means, filted_mu_vals)]

    # band!(ax2, lower_pts, upper_pts, color=(:navyblue, 0.1), transparency=true)
    lines!(ax2, lower_pts, color=:black, linestyle=Linestyle([0, 12, 22]), linewidth=0.9)
    lines!(ax2, upper_pts, color=:black, linestyle=Linestyle([0, 12, 22]), linewidth=0.9)
    scatter!(ax2, mean_pts, color=(:teal, 0.5), marker = :circle)

    ylims!(ax2, 0.6, 3.2)
    zlims!(ax2, 0.0, 0.32)

    colgap!(fig.layout, 10)
    save(joinpath(pic_dir, "nuclear_modes_SetII.pdf"), fig, px_per_unit = 2)
    fig
end

plot_nuclear_modes_transformed()


## figure 3

begin
    function only_steady_state!(ax; target_sets=[], set="SetII", sampling="AA", colors = my_colors)
        parent_dir = pwd()
        all_mus = Set{Float64}()

        for t_set in target_sets
            set_dir = endswith(t_set, "lamda_0.2000") ?
                joinpath(parent_dir, sampling, t_set) :
                joinpath(parent_dir, sampling, t_set, set)
            [push!(all_mus, extract_val("mu", basename(m))) for m in glob("mu*", set_dir) if occursin(r"mu_", m)]
        end

        mu_vals = filter(x -> x >= -2, sort(collect(all_mus)))
        steady_state_curr = fill(NaN, (length(target_sets), length(mu_vals)))

        for (idx_form, t_set) in enumerate(target_sets)
            set_dir = joinpath(parent_dir, sampling, t_set, set)
            for m_dir in glob("mu*", set_dir)
                occursin(r"mu_", m_dir) || continue
                idx_mu = find_nearest(mu_vals, extract_val("mu", basename(m_dir)))
                !isnothing(idx_mu) && (steady_state_curr[idx_form, idx_mu] = get_steady_state(m_dir))
            end
        end

        for i in 1:length(target_sets)
            scatterlines!(ax, mu_vals, steady_state_curr[i,:], color = colors[i])
        end
        return ax
    end

    fig = Figure(size=(800, 350))
    ax1 = Axis(fig[1,1], xlabel=L"\mu_R/\Gamma", ylabel=L"I_{SS}", ylabelrotation=0)
    target_sets=[
        "scan_mu_nuclear_lamda_Qc_01",
        "scan_mu_nuclear_lamda_Qc_05",
        "scan_mu_nuclear_lamda_Qc_15",
        "scan_mu_nuclear_lamda_Qc_20",
    ];
    colors = cgrad(:coolwarm, length(target_sets), categorical = true)
    lam_mins = [parse(Float64, match(r"\d+$", s).match) for s in target_sets] ./ 100.0
    only_steady_state!(ax1; target_sets=target_sets, colors=colors)
    plot_color_bar(ax1, lam_mins; label=L"\lambda_{\mathrm{min}}",
        myticks=[0.0, 0.1, 0.2], loc=(0.8, 0.15))


    ql, qr, lmin, lmax = (-2.5, 2.5, 1.0, 3.0)

    ax_inset = Axis(fig[1, 1],
        width = Relative(0.35), height = Relative(0.35),
        halign = 0.05, valign = 0.95,
        backgroundcolor = :white,
        xticks = [ql, 0.0, qr],
        xgridvisible = false,
        ygridvisible = false,
        rightspinevisible = false,
        topspinevisible = false,
        leftspinevisible = false,
        yticklabelsvisible = false,
        yticksvisible = false,
        xlabelsize = 15,
        )

    x_labels = Dict(ql => L"Q_{c,l}", qr => L"Q_{c,r}", 0.0 => "0")
    ax_inset.xtickformat = values -> [get(x_labels, v, L"") for v in values]

    lines!(ax_inset, [ql-1, ql, 0.0, qr, qr+1], [lmin, lmin, lmax, lmin, lmin],
            color = :blue, linewidth = 2.5)
    lines!(ax_inset, [ql, ql, qr, qr], [0.0, lmin, lmin, 0],
            color=:gray, linestyle=:dash)
    vlines!(ax_inset, 0, color = :black, linewidth = 1.5)

    text!(ax_inset, 0.0, lmax, text = L"\lambda_{\mathrm{max}}",
            align = (:right, :top), offset = (-8, 10), fontsize=16)
    text!(ax_inset, 0.0, lmin, text = L"\lambda_{\mathrm{min}}",
            align = (:right, :bottom), offset = (-2, 0), fontsize=16)
    text!(ax_inset, 0.0, lmax+0.5, text = L"\lambda (Q_c)",
            align = (:center, :top), offset = (20, 10), fontsize=16)
    limits!(ax_inset, -3.6, 3.6, 0, lmax+1)


    ax2 = Axis(fig[1,2], xlabel=L"\mu_R/\Gamma")
    target_sets=[
            "scan_mu_nuclear_lamda_Qt_cos_0",
            "scan_mu_nuclear_lamda_Qt_cos_1",
            "scan_mu_nuclear_lamda_Qt_cos_2",
            "scan_mu_nuclear_lamda_Qt_cos_3",
            "scan_mu_nuclear_lamda_Qt_cos_4",
            "scan_mu_nuclear_lamda_Qt_cos_5",
            "scan_mu_nuclear_lamda_Qt_cos_6",
            "scan_mu_nuclear_lamda_Qt_cos_10",
            ];
    colors = cgrad(:coolwarm, length(target_sets), categorical = true)
    only_steady_state!(ax2; target_sets=target_sets, colors=colors)
    n_vals = [parse(Int, match(r"\d+$", s).match) for s in target_sets]
    plot_color_bar(ax2, n_vals; label=L"n", loc=(0.8, 0.15))

    save(joinpath(pic_dir, "lamda_Q_dependence.pdf"), fig)

    fig
end

## figure 4

begin
    local_size=(1450, 350)
    local_set = "SetII"
    fig = Figure(size=local_size)
    ax11 = Axis(fig[1,1],
        ylabel=L"I_{SS}",
        ylabelrotation=0,
        )
    colors = cgrad(:coolwarm, 9, categorical = true)
    _, lamda_vals = plot_steady_state!(ax11; set=local_set,
        target_dir="scan_mu", left_mu_limit=-3.0, colors=colors)
    text!(ax11, 0.15, 0.0, text="without Nuclear Modes", space=:relative)
    plot_color_bar(ax11, lamda_vals; loc = (0.1, 0.7),
        myticks = [0.0, 0.1, 0.2], width=0.45)
    _, lamda_vals = plot_steady_state!(ax11; set=local_set,
        target_dir=joinpath(pwd(), "juliaHEOM", "scan_mu"),
        left_mu_limit=-3.0, colors=colors, linestyle=:dash)

    ax12 = Axis(fig[1,2],)
    _, lamda_vals = plot_steady_state!(ax12; target_dir="scan_mu_nuclear",
        set=local_set, left_mu_limit=-3.0, colors=colors)

    ax13 = Axis(fig[1,3],)
    plot_steady_state!(ax13; target_dir="scan_mu_nuclear_bath",
        set=local_set, left_mu_limit=-3.0, colors=colors)
    text!(ax13, 0.8, 0.1, text=latexstring("\\eta = \\Delta E/10"),
        space=:relative, align = (:center, :center))

    ax14 = Axis(fig[1,4],)
    colors = cgrad(:coolwarm, 9, categorical = true)
    plot_steady_state!(ax14; target_dir="scan_mu_nuclear_bath2",
        set=local_set, left_mu_limit=-3.0, colors=colors, plot_x="lamda",)
    text!(ax14, 0.8, 0.1, text=latexstring("\\eta = \\Delta E/4"),
        space=:relative, align = (:center, :center))
    # plot_color_bar(ax14, lamda_vals)
    text!(ax14, 0.2, 0.85, text="Set II", space=:relative, font=:bold)

    text!(ax11, 0.015, 0.85, text="(A)", space=:relative, font=:bold)
    text!(ax12, 0.015, 0.85, text="(B)", space=:relative, font=:bold)
    text!(ax13, 0.015, 0.85, text="(C)", space=:relative, font=:bold)
    text!(ax14, 0.015, 0.85, text="(D)", space=:relative, font=:bold)

    linkxaxes!(ax12, ax13, ax14)
    Label(fig[2, :], L"\mu_R/\Gamma", fontsize=30)

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 5)
    save(joinpath(pic_dir, "compare_scan_mu_bath_SetII_mu.pdf"), fig)
    fig
end
