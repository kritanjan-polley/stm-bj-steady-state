# %%
using CairoMakie
using LaTeXStrings
using CSV, DataFrames
using HDF5
using Printf
using Glob
using Statistics

# %%
loadtxt(str) = Matrix(CSV.read(str, DataFrame; delim=' ', header=false, ignorerepeated=true))
find_nearest(array :: AbstractVector, value :: Float64) = argmin(abs.(array .- value))
extract_val(prefix, str) = parse(Float64, match(Regex("$(prefix)_(-?\\d*\\.?\\d+)"), str).captures[1])

pic_dir = joinpath(pwd(), "pic_dir")
mkpath(pic_dir)
my_colors = cgrad(:glasbey_bw_minc_20_maxl_70_n256, 256, categorical=true)
save_fig = true;

# %%

# %%
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

# %%

# %%
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
            p1..., p2..., len/1.3, len/1.3,
            0, false, false,
        )
    end

    path = lift(arc) do arc
        BezierPath([
            MoveTo(points[][1]), arc,
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
            align = (:center, :bottom),
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

        left_bottom_mid_x = X[1] - bar_width / 2
        right_bottom_mid_x = X[2] + bar_width / 2

        x_bar1 = -battery_gap / 2
        x_bar2 = battery_gap / 2

        x_left_stop = x_bar1
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
    plot_x="lamda", plot_every=1, left_mu_limit = -5.0, right_mu_limit=7.0,
    colors = my_colors, if_colorbar=false, linestyle=nothing)

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
            if left_mu_limit <= mu_vals[i] <= right_mu_limit
                k += 1
                valid_mu_len = length(findall(x -> (left_mu_limit <= x <= right_mu_limit), mu_vals))
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
        valid_mu_idx = findall(x -> (left_mu_limit <= x <= right_mu_limit), mu_vals)
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

    return ax1, plot_x == "lamda" ? lamda_vals : mu_vals[findall(x -> (left_mu_limit <= x <= right_mu_limit), mu_vals)]
end

function plot_color_bar(ax, vals;
        loc=(0.8,0.15), myticks=nothing, label=L"\mu_R/\Gamma",
        tickformat=nothing, colors=nothing, width=0.5)

    colors = isnothing(colors) ? cgrad(:coolwarm, length(vals), categorical=true) : colors

    kwargs = (
        bbox=ax.scene.viewport,
        colormap=colors,
        limits=extrema(vals),
        ticklabelsize=18,
        vertical=false,
        halign=loc[1], valign=loc[2],
        width=Relative(width),
        label=label, labelpadding=-70,
    )

    !isnothing(myticks) && (kwargs=(kwargs..., ticks=myticks))
    !isnothing(tickformat) && (kwargs=(kwargs..., tickformat=tickformat))

    Colorbar(fig; kwargs...)
end

# %%

# %%
set_top = "SetI"
set_bottom = "SetIII"

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
        pot1, pot2 = adiabatic_potentials(qt, qc)
        E_minus[i, j] = pot1
        E_plus[i, j]  = pot2
    end
end

top_col_width = 0.8
colors = cgrad(:coolwarm, 9, categorical = true)

fig = Figure(size=(1400, 550))

ax1 = Axis3(fig[1:3, 1],
        xlabel = L"Q_t",
        ylabel = L"Q_c",
        zlabel = "",
        azimuth = deg2rad(45),
        elevation = deg2rad(15),
        aspect=(2,2,2),
        xticks=LinearTicks(3),
        yticks=LinearTicks(3),
        title=L"\textrm{Adiabatic Potential}\, \mathbf{h}(Q_c,Q_t)",
        titlegap=-25,
        )

surface!(ax1, Qt_range, Qc_range, E_minus, colormap=:davos, alpha=0.8, rasterize = 5)
surface!(ax1, Qt_range, Qc_range, E_plus, colormap=:acton, alpha=0.8, rasterize = 5)

# Column 2
g12 = GridLayout(fig[1, 2])
colsize!(g12, 1, Relative(top_col_width))
ax2 = Axis(g12[1, 1])
plot_electrode_diagram!(ax2; voltage=1.8)

ax3 = Axis(fig[2, 2], ylabel=L"I_{SS}",  ylabelrotation=0)
_, lamda_vals = plot_steady_state!(ax3; plot_x ="lamda",
    target_dir="scan_mu",
    set=set_top,
    colors=colors,
    left_mu_limit=-3.0,
    right_mu_limit=5.0,)
_, lamda_vals = plot_steady_state!(ax3; set=set_top, plot_x ="lamda",
    target_dir=joinpath(pwd(), "juliaHEOM", "scan_mu"),
    left_mu_limit=-3.0, right_mu_limit=5.0,
    colors=colors, linestyle=:dash,)

ax8 = Axis(fig[3, 2], ylabel=L"I_{SS}", ylabelrotation=0)
_, lamda_vals = plot_steady_state!(ax8; plot_x ="lamda",
    target_dir="scan_mu",
    set=set_bottom,
    colors=colors,
    left_mu_limit=-3.0,
    right_mu_limit=5.0,
    )
plot_color_bar(ax8, lamda_vals; loc=(0.1, 0.8), label=L"\lambda/\Gamma")
_, lamda_vals = plot_steady_state!(ax8; set=set_bottom, plot_x ="lamda",
    target_dir=joinpath(pwd(), "juliaHEOM", "scan_mu"),
    left_mu_limit=-3.0, right_mu_limit=5.0,
    colors=colors, linestyle=:dash,)

# Column 3
g13 = GridLayout(fig[1, 3])
colsize!(g13, 1, Relative(top_col_width))
ax4 = Axis(g13[1, 1])
plot_electrode_diagram!(ax4; voltage=1.4, if_oscillator=true)

ax5 = Axis(fig[2, 3])
_, lamda_vals = plot_steady_state!(ax5; plot_x ="lamda",
    set=set_top,
    colors=colors,
    left_mu_limit=-3.0,
    target_dir="scan_mu_nuclear_temp",
    )


ax9 = Axis(fig[3, 3])
_, lamda_vals = plot_steady_state!(ax9; plot_x ="lamda",
    set=set_bottom,
    colors=colors,
    left_mu_limit=-3.0,
    target_dir="scan_mu_nuclear_temp",
)

# Column 4
g14 = GridLayout(fig[1, 4])
colsize!(g14, 1, Relative(top_col_width))
ax6 = Axis(g14[1, 1])
plot_electrode_diagram!(ax6; voltage=0.8, if_oscillator=true, connect_both=true)

ax7 = Axis(fig[2, 4])
_, lamda_vals = plot_steady_state!(ax7; plot_x ="lamda",
    target_dir="scan_mu_nuclear_check_coupled_temp",
    set=set_top,
    plot_every=1,
    colors=colors[1:2:end],)
text!(ax7, 0.6, 0.2; text="Set I", space=:relative, font=:bold)

ax10 = Axis(fig[3, 4])
_, lamda_vals = plot_steady_state!(ax10; plot_x ="lamda",
    target_dir="scan_mu_nuclear_check_coupled_temp",
    set=set_bottom,
    plot_every=1,
    colors=colors[1:2:end])
text!(ax10, 0.6, 0.2; text="Set III", space=:relative, font=:bold)

linkxaxes!(ax3, ax8)
linkxaxes!(ax5, ax9)
linkxaxes!(ax7, ax10)

hidexdecorations!(ax3; grid=false, ticks=false)
hidexdecorations!(ax5; grid=false, ticks=false)
hidexdecorations!(ax7; grid=false, ticks=false)

Label(fig[1, 1, TopLeft()], "(A)", font=:bold, padding=(0, 20, -50, 0), halign=:right)
Label(fig[1, 2, TopLeft()], "(B)", font=:bold, padding=(0, 5, -50, 0), halign=:right)
Label(fig[1, 3, TopLeft()], "(C)", font=:bold, padding=(0, 5, -50, 0), halign=:right)
Label(fig[1, 4, TopLeft()], "(D)", font=:bold, padding=(0, 5, -50, 0), halign=:right)

Label(fig[4, 2:4], L"\mu_R/\Gamma", tellwidth=false, fontsize=28)

rowsize!(fig.layout, 1, Auto(0.45))
rowsize!(fig.layout, 2, Auto(0.8))
rowsize!(fig.layout, 3, Auto(0.8))

rowgap!(fig.layout, 1, 5)
rowgap!(fig.layout, 2, 5)
colgap!(fig.layout, 1, 0)

if save_fig
    filename = "nuclear_pot_system_compare_$(set_top)_$(set_bottom).pdf"
    @info "Saving $(filename)"
    save(joinpath(pic_dir, filename), fig)
end

fig

# %%

# %%
function plot_auto!(ax; kwargs...)
    _, mu = plot_steady_state!(ax; colors=my_colors, kwargs...)
    empty!(ax)
    cols = cgrad(:coolwarm, length(mu), categorical=true)
    plot_steady_state!(ax; colors=cols, kwargs...)
    mu, cols
end

begin
    println("generating Figure 1")

    set_top, set_bottom = "SetI", "SetIII"
    delta = 1.0
    omega_t, omega_c = 0.1delta, 0.2delta
    epsilon1, epsilon2 = -delta/2, delta/2
    kappa1, kappa2, gamma = -0.2delta, 0.3delta, 0.05delta

    function adiabatic_potentials(Qt,Qc)
        V0 = 0.5omega_t*Qt^2 + 0.5omega_c*Qc^2
        V11, V22, V12 = V0+epsilon1+kappa1*Qt, V0+epsilon2+kappa2*Qt, gamma*Qc
        m, d = 0.5(V11+V22), 0.5(V11-V22)
        g = sqrt(d^2+V12^2)
        m-g, m+g
    end

    nr = 101
    Qt, Qc = range(-4,0,length=nr), range(-2,2,length=nr)
    Em, Ep = zeros(nr,nr), zeros(nr,nr)
    for (i,qt) in enumerate(Qt), (j,qc) in enumerate(Qc)
        Em[i,j], Ep[i,j] = adiabatic_potentials(qt,qc)
    end

    fig = Figure(size=(1400,550))
    w = 0.8

    ax1 = Axis3(fig[1:3,1], xlabel=L"Q_t", ylabel=L"Q_c", zlabel="",
        azimuth=deg2rad(45), elevation=deg2rad(15), aspect=(2,2,2),
        xticks=LinearTicks(3), yticks=LinearTicks(3),
        title=L"\textrm{Adiabatic Potential}\,\mathbf{h}(Q_c,Q_t)", titlegap=-25)
    surface!(ax1,Qt,Qc,Em,colormap=:davos,alpha=0.8,rasterize=5)
    surface!(ax1,Qt,Qc,Ep,colormap=:acton,alpha=0.8,rasterize=5)

    g12=GridLayout(fig[1,2]); colsize!(g12,1,Relative(w))
    g13=GridLayout(fig[1,3]); colsize!(g13,1,Relative(w))
    g14=GridLayout(fig[1,4]); colsize!(g14,1,Relative(w))

    ax2=Axis(g12[1,1]); plot_electrode_diagram!(ax2; voltage=1.8)
    ax4=Axis(g13[1,1]); plot_electrode_diagram!(ax4; voltage=1.4, if_oscillator=true)
    ax6=Axis(g14[1,1]); plot_electrode_diagram!(ax6; voltage=0.8, if_oscillator=true, connect_both=true)

    ax3=Axis(fig[2,2],ylabel=L"I_{SS}",ylabelrotation=0)
    ax8=Axis(fig[3,2],ylabel=L"I_{SS}",ylabelrotation=0)
    ax5=Axis(fig[2,3]); ax9=Axis(fig[3,3])
    ax7=Axis(fig[2,4]); ax10=Axis(fig[3,4])

    # Probe all datasets for their actual mu values
    function get_mu!(ax; kwargs...)
        _, mu = plot_steady_state!(ax; colors=my_colors, kwargs...)
        empty!(ax)
        mu
    end

    mu3 = get_mu!(ax3;  plot_x="mu",set=set_top,target_dir="scan_mu",left_mu_limit=-3.0,right_mu_limit=5.0)
    mu8 = get_mu!(ax8;  plot_x="mu",set=set_bottom,target_dir="scan_mu",left_mu_limit=-3.0,right_mu_limit=5.0)
    mu5 = get_mu!(ax5;  plot_x="mu",set=set_top,target_dir="scan_mu_nuclear_temp",left_mu_limit=-3.0,right_mu_limit=5.0)
    mu9 = get_mu!(ax9;  plot_x="mu",set=set_bottom,target_dir="scan_mu_nuclear_temp",left_mu_limit=-3.0,right_mu_limit=5.0)
    mu7 = get_mu!(ax7;  plot_x="mu",set=set_top,target_dir="scan_mu_nuclear_check_coupled_temp",plot_every=1)
    mu10 = get_mu!(ax10; plot_x="mu",set=set_bottom,target_dir="scan_mu_nuclear_check_coupled_temp",plot_every=1)

    mu_ref = sort(unique(vcat(mu3,mu8,mu5,mu9,mu7,mu10)))
    cmap, mumin, mumax = cgrad(:coolwarm), minimum(mu_ref), maximum(mu_ref)
    getcolors(mu) = [cmap[clamp((x-mumin)/(mumax-mumin),0,1)] for x in mu]

    c3,c8,c5,c9,c7,c10 = getcolors(mu3),getcolors(mu8),getcolors(mu5),
                          getcolors(mu9),getcolors(mu7),getcolors(mu10)

    # Column B
    plot_steady_state!(ax3; plot_x="mu",set=set_top,target_dir="scan_mu",
        colors=c3,left_mu_limit=-3.0,right_mu_limit=5.0)
    plot_steady_state!(ax3; plot_x="mu",set=set_top,
        target_dir=joinpath(pwd(),"juliaHEOM","scan_mu"),colors=c3,
        left_mu_limit=-3.0,right_mu_limit=5.0,linestyle=:dash)

    plot_steady_state!(ax8; plot_x="mu",set=set_bottom,target_dir="scan_mu",
        colors=c8,left_mu_limit=-3.0,right_mu_limit=5.0)
    plot_steady_state!(ax8; plot_x="mu",set=set_bottom,
        target_dir=joinpath(pwd(),"juliaHEOM","scan_mu"),colors=c8,
        left_mu_limit=-3.0,right_mu_limit=5.0,linestyle=:dash)

    # Column C
    plot_steady_state!(ax5; plot_x="mu",set=set_top,target_dir="scan_mu_nuclear_temp",
        colors=c5,left_mu_limit=-3.0,right_mu_limit=5.0)
    plot_steady_state!(ax9; plot_x="mu",set=set_bottom,target_dir="scan_mu_nuclear_temp",
        colors=c9,left_mu_limit=-3.0,right_mu_limit=5.0)

    # Column D
    plot_steady_state!(ax7; plot_x="mu",set=set_top,
        target_dir="scan_mu_nuclear_check_coupled_temp",colors=c7,plot_every=1)
    plot_steady_state!(ax10; plot_x="mu",set=set_bottom,
        target_dir="scan_mu_nuclear_check_coupled_temp",colors=c10,plot_every=1)

    text!(ax7,0.6,0.2,text="Set I",space=:relative,font=:bold)
    text!(ax10,0.6,0.2,text="Set III",space=:relative,font=:bold)

    # for ax in (ax3, ax8, ax5, ax9, ax7, ax10)
        plot_color_bar(ax8, mu_ref; loc=(0.1,0.8), label=L"\mu_R/\Gamma",
        myticks= range(minimum(mu_ref), maximum(mu_ref), length=3),
        tickformat=x -> [@sprintf("%.1f", val) for val in x])
    # end

    linkxaxes!(ax3,ax8); linkxaxes!(ax5,ax9); linkxaxes!(ax7,ax10)
    for ax in (ax3,ax5,ax7)
        hidexdecorations!(ax; grid=false,ticks=false)
    end

    for (i,p) in enumerate(((0,20,-50,0),(0,5,-50,0),(0,5,-50,0),(0,5,-50,0)))
        Label(fig[1,i,TopLeft()],"($(Char(64+i)))",font=:bold,padding=p,halign=:right)
    end

    Label(fig[4,2:4],L"\lambda/\Gamma",tellwidth=false,fontsize=28)
    rowsize!(fig.layout,1,Auto(0.45)); rowsize!(fig.layout,2,Auto(0.8)); rowsize!(fig.layout,3,Auto(0.8))
    rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,5); colgap!(fig.layout,1,0)

    # filename = "current_vs_lambda_$(set_top)_$(set_bottom).pdf"
    # @info "Saving $filename"
    # save(joinpath(pic_dir,filename),fig)
    fig
end

# %%

# %%

# %%

# %%
function plot_nuclear_modes_transformed!(
    figpos;
    sampling="AA",
    target_dir="scan_mu_nuclear",
    set="SetII",
    lamda_dir="lamda_0.0750",
    mu_min=0.0,
    mu_max=3.0,
    z_max = 0.6,
    title="",
    )
    data_dir = joinpath(pwd(), sampling, target_dir, set, lamda_dir)
    mu_dirs = glob("mu_*", data_dir)

    pairs = [
        (parse(Float64, match(r"mu_(-?\d+\.?\d*)", basename(p)).captures[1]), p)
        for p in mu_dirs
    ]

    pairs = sort(filter(((mu, p),) -> mu_min <= mu <= mu_max, pairs), by=first)

    filtered_mu_vals = first.(pairs)
    filtered_mu_dirs = last.(pairs)

    raw_qc_data = Dict{Int,Vector{Float64}}()
    raw_qt_data = Dict{Int,Vector{Float64}}()

    for (i, mu) in enumerate(filtered_mu_vals)
        datafile = joinpath(filtered_mu_dirs[i], "qc_qt_average.h5")
        @assert isfile(datafile)

        h5open(datafile, "r") do f
            raw = read(f["data"])
            raw_qc_data[i] = raw[:, 1]
            raw_qt_data[i] = raw[:, 2]
        end
    end

    ax1 = Axis3(figpos[1, 1],
        xlabel=L"Q_c",
        ylabel=L"\mu_R/\Gamma",
        zlabel="",
        yticks=LinearTicks(5),
        zticks=LinearTicks(4),
        azimuth=deg2rad(-49),
        elevation=deg2rad(25),
        aspect=(3, 3, 1.5),
        xgridcolor=(:black, 0.05),
        ygridcolor=(:black, 0.05),
    )

    for (i, mu) in enumerate(filtered_mu_vals)
        plt = density!(
            ax1,
            raw_qc_data[i],
            strokewidth=1.2,
            color=(:orchid4, 0.2),
            strokecolor=:orchid4,
        )
        rotate!(plt, Vec3f(1, 0, 0), pi / 2)
        translate!(plt, 0, mu, 0)
    end

    text!(
        ax1,
        0.25, 0.75,
        text=title,
        space=:relative,
        align=(:left, :top),
        fontsize=24,
        font=:bold,
    )

    ylims!(ax1, mu_min - 0.2, mu_max + 0.2)
    zlims!(ax1, 0.0, z_max)

    ax2 = Axis3(figpos[1, 2],
        xlabel=L"Q_t",
        ylabel=L"\mu_R/\Gamma",
        zlabel="",
        yticks=LinearTicks(5),
        zticks=LinearTicks(4),
        azimuth=deg2rad(49),
        elevation=deg2rad(25),
        aspect=(3, 3, 1.5),
        xgridcolor=(:black, 0.05),
        ygridcolor=(:black, 0.05),
    )

    for (i, mu) in enumerate(filtered_mu_vals)
        plt = density!(
            ax2, raw_qt_data[i],
            strokewidth=1.2,
            color=(:turquoise4, 0.2),
            strokecolor=:turquoise4,
        )
        rotate!(plt, Vec3f(1, 0, 0), pi / 2)
        translate!(plt, 0, mu, 0)
    end


    qt_means = [mean(raw_qt_data[i]) for i in eachindex(filtered_mu_vals)]
    qt_stds  = [std(raw_qt_data[i]) for i in eachindex(filtered_mu_vals)]

    z_floor = 0.001

    lower_pts = [
        Point3f(m - s, y, z_floor)
        for (m, s, y) in zip(qt_means, qt_stds, filtered_mu_vals)
    ]

    upper_pts = [
        Point3f(m + s, y, z_floor)
        for (m, s, y) in zip(qt_means, qt_stds, filtered_mu_vals)
    ]

    mean_pts = [
        Point3f(m, y, z_floor)
        for (m, y) in zip(qt_means, filtered_mu_vals)
    ]

    lines!(
        ax2, [-2.0, -2.0], [mu_min - 0.2, mu_max + 0.2], [0.0, 0.0],
        color=:black, linestyle=:dash, linewidth=1.0,
    )

    lines!(ax2, lower_pts, color=:black, linestyle=Linestyle([0, 12, 22]), linewidth=0.9)
    lines!(ax2, upper_pts, color=:black, linestyle=Linestyle([0, 12, 22]), linewidth=0.9)
    scatter!(ax2, mean_pts, color=(:teal, 0.5), marker=:circle, markersize=15)

    ylims!(ax2, mu_min - 0.2, mu_max + 0.2)
    zlims!(ax2, 0.0, z_max/2)

    return ax1, ax2
end

function plot_nuclear_modes_SetI_SetIII()
    fig = Figure(size=(1050, 655), figure_padding=(0, -10, 0, 0))

    plot_nuclear_modes_transformed!(
        fig[1, 1];
        target_dir="scan_mu_nuclear_temp",
        set="SetI",
        lamda_dir="lamda_0.1000",
        title="Set I", z_max = 0.65,
        mu_min = -3.0, mu_max = 5.0,
    )

    plot_nuclear_modes_transformed!(
        fig[2, 1];
        target_dir="scan_mu_nuclear_temp",
        set="SetIII",
        lamda_dir="lamda_0.1000",
        title="Set III", z_max = 0.35,
        mu_min = -3.0, mu_max = 5.0,
    )

    rowgap!(fig.layout, -50)
    colgap!(fig.layout, 10)

    if save_fig
        filename = "nuclear_modes_SetI_SetIII.pdf"
        @info "Saving $(filename)"
        save(joinpath(pic_dir, filename), fig, px_per_unit=2)
    end

    fig
end

plot_nuclear_modes_SetI_SetIII()

# %%

# %%
function only_steady_state!(ax; target_sets=[], set="SetII", sampling="AA", colors=my_colors)
    parent_dir = pwd()

    for (idx_form, t_set) in enumerate(target_sets)
        set_dir = endswith(t_set, "lamda_0.2000") ?
            joinpath(parent_dir, sampling, t_set) :
            joinpath(parent_dir, sampling, t_set, set)

        mus = Float64[]
        steady_state_curr = Float64[]

        for m_dir in glob("mu*", set_dir)
            occursin(r"mu_", m_dir) || continue

            mu_val = extract_val("mu", basename(m_dir))
            mu_val >= -2.0 || continue

            try
                curr = get_steady_state(m_dir)
                push!(mus, mu_val)
                push!(steady_state_curr, curr)
            catch err
                @warn "Skipping failed directory" m_dir err
            end
        end

        perm = sortperm(mus)

        scatterlines!(
            ax, mus[perm], steady_state_curr[perm], color=colors[idx_form],
            )
    end
    return ax
end

local_set = "SetII"
fig = Figure(size=(800, 350))
ax1 = Axis(fig[1,1], xlabel=L"\mu_R/\Gamma", ylabel=L"I_{SS}", ylabelrotation=0)
target_sets=[
    "scan_mu_nuclear_lamda_Qc_01_temp",
    "scan_mu_nuclear_lamda_Qc_05_temp",
    "scan_mu_nuclear_lamda_Qc_15_temp",
    "scan_mu_nuclear_lamda_Qc_20_temp",
];
lam_mins = [parse(Float64, match(r"Qc_(\d+)_", s).captures[1]) for s in target_sets] ./ 100.0
colors = cgrad(:coolwarm, length(target_sets), categorical = true)
only_steady_state!(ax1; target_sets=target_sets, set=local_set, colors=colors)
plot_color_bar(ax1, lam_mins; label=L"\lambda_{\mathrm{min}}",
    myticks=[0.0, 0.1, 0.2], loc=(0.8, 0.15))


ql, qr, lmin, lmax = (-2.5, 2.5, 1.0, 3.0)

ax_inset = Axis(fig[1, 1],
    width = Relative(0.35), height = Relative(0.35),
    halign = 0.05, valign = 0.95,
    backgroundcolor = :white,
    xticks = [ql, 0.0, qr],
    xgridvisible = false, ygridvisible = false,
    rightspinevisible = false,
    topspinevisible = false,
    leftspinevisible = false,
    yticklabelsvisible = false, yticksvisible = false,
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
        "scan_mu_nuclear_lamda_Qt_cos_1_temp",
        "scan_mu_nuclear_lamda_Qt_cos_2_temp",
        "scan_mu_nuclear_lamda_Qt_cos_3_temp",
        "scan_mu_nuclear_lamda_Qt_cos_4_temp",
        "scan_mu_nuclear_lamda_Qt_cos_5_temp",
        "scan_mu_nuclear_lamda_Qt_cos_6_temp",
        "scan_mu_nuclear_lamda_Qt_cos_10_temp",
        ];
colors = cgrad(:coolwarm, length(target_sets), categorical = true)
only_steady_state!(ax2; target_sets=target_sets, colors=colors)
n_vals = [parse(Int, match(r"cos_(\d+)", s).captures[1]) for s in target_sets]
plot_color_bar(ax2, n_vals; label=L"n", loc=(0.8, 0.15), myticks=[1, 5, 10],)
text!(ax2, 0.45, 1.0, text="Set II",
    space=:relative, align=(:left, :top),
    fontsize=24, font=:bold,
)

ax_cos = Axis(fig[1,2],
    width=Relative(0.38), height=Relative(0.3),
    halign=0.05, valign=0.95,
    backgroundcolor=:white,
    xticks=[-pi, 0, pi],
    xtickformat = _ -> [L"Q_{c,l}", L"0", L"Q_{c,r}"],
    yticksvisible=false, yticklabelsvisible=false,
    xgridvisible=false, ygridvisible=false,
    rightspinevisible=false, topspinevisible=false,
    leftspinevisible=false,
    xticklabelsize=18,
)

x = range(-pi, pi, length=300)

for n in 1:2
    lines!(ax_cos, x, cos.(n .* x),
        color=colors[n], linewidth=2)
end

vlines!(ax_cos, 0.0, color=:black, linewidth=1.5)
text!(ax_cos, 1.5, 1.2, text=L"\lambda(Q_t)",
    align=(:center, :top), fontsize=18, font=:bold)

if save_fig
    filename = "lamda_Q_dependence_$(local_set).pdf"
    @info "Saving $(filename)"
    save(joinpath(pic_dir, filename), fig)
end

fig

# %%

# %%
local_size=(1450, 350)
local_set = "SetII"
fig = Figure(size=local_size)
colors = cgrad(:coolwarm, 9, categorical = true)

ax12 = Axis(fig[1,1], ylabel=L"I_{SS}", ylabelrotation=0,)
_, lamda_vals = plot_steady_state!(ax12; target_dir="scan_mu_nuclear_temp",
    set=local_set, left_mu_limit=-3.0, colors=colors)
text!(ax12, 0.8, 0.1, text=latexstring("\\eta = 0"),
    space=:relative, align = (:center, :center))
plot_color_bar(ax12, lamda_vals; loc = (0.1, 0.7),
    myticks = [0.0, 0.1, 0.2], width=0.45, label=L"\lambda/\Gamma")

ax13 = Axis(fig[1,2],)
plot_steady_state!(ax13; target_dir="scan_mu_nuclear_bath_quarter_temp",
    set=local_set, left_mu_limit=-3.0, colors=colors)
text!(ax13, 0.8, 0.1, text=latexstring("\\eta = \\Delta/4"),
    space=:relative, align = (:center, :center))

ax14 = Axis(fig[1,3],)
plot_steady_state!(ax14; target_dir="scan_mu_nuclear_bath_one_temp",
    set=local_set, left_mu_limit=-3.0, colors=colors, plot_x="lamda",)
text!(ax14, 0.8, 0.1, text=latexstring("\\eta = \\Delta"),
    space=:relative, align = (:center, :center))

ax15 = Axis(fig[1,4],)
plot_steady_state!(ax15; target_dir="scan_mu_nuclear_bath_two_temp",
    set=local_set, left_mu_limit=-3.0, colors=colors, plot_x="lamda",)
text!(ax15, 0.8, 0.1, text=latexstring("\\eta = 2\\Delta"),
    space=:relative, align = (:center, :center))
text!(ax15, 0.1, 0.7, text="Set II", space=:relative, font=:bold)

text!(ax12, 0.015, 0.85, text="(A)", space=:relative, font=:bold)
text!(ax13, 0.015, 0.85, text="(B)", space=:relative, font=:bold)
text!(ax14, 0.015, 0.85, text="(C)", space=:relative, font=:bold)
text!(ax15, 0.015, 0.85, text="(D)", space=:relative, font=:bold)

linkxaxes!(ax12, ax13, ax14, ax15)
Label(fig[2, :], L"\mu_R/\Gamma", fontsize=30)

colgap!(fig.layout, 10)
rowgap!(fig.layout, 5)

if save_fig
    filename = "compare_scan_mu_bath_$(local_set)_mu.pdf"
    @info "Saving $(filename)"
    save(joinpath(pic_dir, filename), fig)
end

fig

# %%

# %%
begin
    local_size = (1450, 350)
    local_set = "SetII"
    fig = Figure(size=local_size)

    function get_mu!(ax; kwargs...)
        _, mu = plot_steady_state!(ax; colors=my_colors, kwargs...)
        empty!(ax)
        mu
    end

    ax12 = Axis(fig[1,1], ylabel=L"I_{SS}", ylabelrotation=0)
    ax13 = Axis(fig[1,2])
    ax14 = Axis(fig[1,3])
    ax15 = Axis(fig[1,4])

    mu12 = get_mu!(ax12; target_dir="scan_mu_nuclear_temp",
        set=local_set, left_mu_limit=-3.0, plot_x="mu")
    mu13 = get_mu!(ax13; target_dir="scan_mu_nuclear_bath_quarter_temp",
        set=local_set, left_mu_limit=-3.0, plot_x="mu")
    mu14 = get_mu!(ax14; target_dir="scan_mu_nuclear_bath_one_temp",
        set=local_set, left_mu_limit=-3.0, plot_x="mu")
    mu15 = get_mu!(ax15; target_dir="scan_mu_nuclear_bath_two_temp",
        set=local_set, left_mu_limit=-3.0, plot_x="mu")

    mu_ref = sort(unique(vcat(mu12, mu13, mu14, mu15)))
    cmap = cgrad(:coolwarm)
    mu0, mu1 = minimum(mu_ref), maximum(mu_ref)
    getcols(mu) = [cmap[clamp((x-mu0)/(mu1-mu0), 0, 1)] for x in mu]

    c12, c13, c14, c15 = getcols(mu12), getcols(mu13), getcols(mu14), getcols(mu15)

    plot_steady_state!(ax12; target_dir="scan_mu_nuclear_temp",
        set=local_set, left_mu_limit=-3.0, colors=c12, plot_x="mu")
    plot_steady_state!(ax13; target_dir="scan_mu_nuclear_bath_quarter_temp",
        set=local_set, left_mu_limit=-3.0, colors=c13, plot_x="mu")
    plot_steady_state!(ax14; target_dir="scan_mu_nuclear_bath_one_temp",
        set=local_set, left_mu_limit=-3.0, colors=c14, plot_x="mu")
    plot_steady_state!(ax15; target_dir="scan_mu_nuclear_bath_two_temp",
        set=local_set, left_mu_limit=-3.0, colors=c15, plot_x="mu")

    text!(ax12, 0.8, 0.5, text=latexstring("\\eta = 0"),
        space=:relative, align=(:center, :center))
    text!(ax13, 0.8, 0.5, text=latexstring("\\eta = \\Delta/4"),
        space=:relative, align=(:center, :center))
    text!(ax14, 0.8, 0.5, text=latexstring("\\eta = \\Delta"),
        space=:relative, align=(:center, :center))
    text!(ax15, 0.8, 0.5, text=latexstring("\\eta = 2\\Delta"),
        space=:relative, align=(:center, :center))
    text!(ax15, 0.1, 0.7, text="Set II", space=:relative, font=:bold)

    for (ax, lab) in zip((ax12, ax13, ax14, ax15), ("(A)", "(B)", "(C)", "(D)"))
        text!(ax, 0.015, 0.85, text=lab, space=:relative, font=:bold)
    end

    cb_ticks = range(mu0, mu1, length=3)
    plot_color_bar(ax12, mu_ref; loc=(0.1, 0.7), width=0.45,
        label=L"\mu_R/\Gamma", colors=getcols(mu_ref),
        myticks=cb_ticks, tickformat=x -> [@sprintf("%.1f", v) for v in x])

    linkxaxes!(ax12, ax13, ax14, ax15)
    Label(fig[2, :], L"\lambda/\Gamma", fontsize=30)

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 5)

    if save_fig
        filename = "compare_scan_mu_bath_$(local_set)_lambda.pdf"
        @info "Saving $(filename)"
        save(joinpath(pic_dir, filename), fig)
    end

    fig
end
