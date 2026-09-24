using CairoMakie
using LaTeXStrings
using CSV, DataFrames
using Printf
using Glob
using Statistics

loadtxt(str::AbstractString) = Matrix(CSV.read(
        str, DataFrame; delim=' ', header=false, ignorerepeated=true
        ))
find_nearest(array :: AbstractVector, value :: Float64) = argmin(abs.(array .- value))
extract_val(prefix, str) = parse(Float64,
            match(Regex("$(prefix)_(-?\\d*\\.?\\d+)"), str).captures[1])

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
    palette = (color = my_colors,),
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


function plot_steady_state_only!(ax; set="SetII", target_dir="scan_mu",
                                 plot_x="lamda", plot_every=1, left_mu_limit = -5.0,
                                 colors = my_colors, if_colorbar=false)

    target_dir = isabspath(target_dir) ? target_dir : joinpath(pwd(), "AA", target_dir)
    set_dir = joinpath(target_dir, set)

    lamda_dirs = sort(glob("lamda*", set_dir), by = x -> extract_val("lamda", x))
    lamda_vals = zeros(length(lamda_dirs))

    mu_vals = zeros(length(glob("mu*", lamda_dirs[1])))
    steady_state_curr = zeros((length(lamda_dirs), length(mu_vals)))

    for (l, lamda) in enumerate(lamda_dirs)
        mu_dirs = sort(glob("mu*", lamda), by = x -> extract_val("mu", x))
        lamda_vals[l] = extract_val("lamda", lamda)

        for (m, mu) in enumerate(mu_dirs)
            mu_val = extract_val("mu", mu)
            mu_vals[m] = contains(lowercase(mu), "heom") ? -mu_val : mu_val

            if contains(lowercase(mu), "heom")
                f1 = loadtxt(joinpath(mu, "current.HEOM.txt"))
                current = -0.5 .* (f1[:,2] .- f1[:,3])
            elseif contains(lowercase(mu), "master")
                f1 = loadtxt(joinpath(mu, "current.masterEq.txt"))
                if size(f1, 2) == 2
                    current = f1[:,2]
                else
                    current = 0.5 .* (f1[:,2] .- f1[:,3])
                end
            else
                f1 = loadtxt(joinpath(mu, "current_average.txt"))
                current = -0.5 .* (f1[:,2] .- f1[:,3])
            end

            times = f1[:,1]
            start_idx = find_nearest(times, times[end] - 3.0)
            steady_state_curr[l, m] = mean(current[start_idx:end])
        end
    end

    plt = image!(ax, mu_vals[1]..mu_vals[end], lamda_vals[1]..lamda_vals[end],
                 transpose(steady_state_curr), colormap=:nipy_spectral)
    xlims!(ax, -3.0, 4.0)
    return ax, plt
end

function check_population!(ax; set="SetI", target_dir="scan_mu_nuclear",
                           muL=4.0, lamda=0.2, filetype="current_average.txt")

    base_path(method) = joinpath(pwd(), method, target_dir, set,
                                 @sprintf("lamda_%.4f", lamda),
                                 @sprintf("mu_%.1f", muL), filetype)

    filename_AA, filename_Wi = base_path("AA"), base_path("Wigner")
    @assert isfile(filename_AA) && isfile(filename_Wi)

    f1, f2 = loadtxt(filename_AA), loadtxt(filename_Wi)

    if filetype == "population_average.txt"
        lines!(ax, f1[:,1], f1[:,2], label=L"\rho_{11}")
        lines!(ax, f1[:,1], f1[:,3], label=L"\rho_{22}")
        lines!(ax, f2[:,1], f2[:,2], label=L"\rho_{11}", linestyle=:dash)
        lines!(ax, f2[:,1], f2[:,3], label=L"\rho_{22}", linestyle=:dash)
    elseif filetype == "current_average.txt"
        lines!(ax, f1[:,1], -f1[:,2] .+ f1[:,3],
            label=L"I(t)\, (\mathrm{AA})")
        lines!(ax, f2[:,1], -f2[:,2] .+ f2[:,3],
            label=L"I(t)\, (\mathrm{Wigner})",
            linestyle=:dash)
    elseif filetype == "photon_average.txt"
        lines!(ax, f1[:,1], f1[:,2],
            label=L"\langle a_{\mathrm{cav}}^{\dagger}a_{\mathrm{cav}}\rangle (\mathrm{AA})")
        lines!(ax, f2[:,1], f2[:,2],
            label=L"\langle a_{\mathrm{cav}}^{\dagger}a_{\mathrm{cav}}\rangle (\mathrm{Wigner})",
            linestyle=:dash)
    else
        @error "Unknown filetype"
    end
    return ax
end

function plot_dynamics!(pos, set, muL, heom_file, aa_file; label_str=nothing)
    spectra_set = joinpath(pwd(), "juliaHEOM", "scan_mu", set)
    lamda_dirs = sort(glob("lamda*", spectra_set), by = x -> extract_val("lamda", x))
    heom_mu = -muL

    ax = Axis(pos, xlabel = L"\Gamma t", xlabelsize = 18, ylabelsize = 18)
    lamda_vals = Float64[]

    for (idx, lamda_path) in enumerate(lamda_dirs)
        lamda_val = extract_val("lamda", lamda_path)
        push!(lamda_vals, lamda_val)
        lamda_name = basename(lamda_path)

        # HEOM
        f1 = loadtxt(joinpath(lamda_path, @sprintf("mu_%.1f", heom_mu), heom_file))
        lines!(ax, f1[:,1], f1[:,2], color=local_colors[idx],
            label=latexstring("\\lambda/\\Gamma = $(@sprintf("%.3f", lamda_val))"))

        # semiclassical
        f2 = loadtxt(joinpath(pwd(), "AA", "scan_mu", set,
            lamda_name, @sprintf("mu_%.1f", muL), aa_file))
        lines!(ax, f2[:,1], f2[:,2], color=local_colors[idx], linestyle=:dash)
    end

    if !isnothing(label_str)
        text!(ax, 0.05, 0.3, text = LaTeXString(label_str),
            space=:relative, fontsize=18)
    end
    return ax, lamda_vals
end


## Fig S1

fig = Figure(size=(1000, 300), linewidth=2.5)

configs_s1 = [
    (set="SetIII", mu=-1.0, lamda=0.2, ftype="population_average.txt",
        yl="Population", title_pos=(0.5, 0.8), lab_pos=(0.3, 0.7),
        leg_pos=:rc, banks=2, xmax=nothing),
    (set="SetI",   mu=4.0,  lamda=0.2, ftype="current_average.txt",
        yl="Current",  title_pos=(0.5, 0.5), lab_pos=(0.3, 0.4),
        leg_pos=:rt, banks=1, xmax=nothing),
    (set="SetII",  mu=2.0,  lamda=0.1, ftype="photon_average.txt",
        yl="Photon Number", title_pos=(0.2, 0.9), lab_pos=(0.1, 0.8),
        leg_pos=:rb, banks=1, xmax=20.0)
]

for (i, c) in enumerate(configs_s1)
    ax = Axis(fig[1, i],
        xlabel=L"\Gamma t",
        ylabel=c.yl,
        xtickformat = x -> [@sprintf("%.0f", val) for val in x])
    check_population!(ax; set=c.set, muL=c.mu, lamda=c.lamda, filetype=c.ftype)

    text!(ax, c.title_pos...,
        text="Set $(split(c.set, "Set")[2])",
        space=:relative, fontsize=20)
    text!(ax, c.lab_pos...,
        text=latexstring(
            @sprintf("\\mu_R=%.1f\\Gamma \\,\\lambda=%.1f\\Gamma", c.mu, c.lamda)
            ),
        space=:relative, fontsize=20)

    axislegend(ax, nbanks=c.banks,
        orientation=(c.banks > 1 ? :horizontal : :vertical),
        position=c.leg_pos)
    isnothing(c.xmax) || xlims!(ax, nothing, c.xmax)
end

# save(joinpath(pic_dir, "compare_AA_Wigner.pdf"), fig)
fig

## fig S2

fig = Figure(size=(1200, 350))

configs_s2 = [
    (title="Semiclassical", dir="scan_mu"),
    (title="HEOM", dir=joinpath(pwd(), "juliaHEOM", "scan_mu")),
    # (title="Lindblad", dir=joinpath(pwd(), "masterEq", "scan_mu"),
    (title="Lindblad", dir=joinpath(pwd(), "masterEq2", "Linblad", "bare2LS")),
    # (title="Redfield", dir=joinpath(pwd(), "masterEq2", "Redfield", "bare2LS")),
]

for (i, c) in enumerate(configs_s2)
    ax = Axis(fig[1, 2i-1],
        xlabel = L"\mu_R/\Gamma",
        xticks = LinearTicks(4), title=c.title)
    if i == 1
        ax.ylabel = L"\lambda /\Gamma"
    end

    kwargs = isnothing(c.dir) ? () : (target_dir=c.dir,)
    _, plt = plot_steady_state_only!(ax; set="SetII", kwargs...)

    if i == 1
        text!(ax, 0.2, 0.1, text="Set II", space=:relative, color="white")
    end
    # xlims!(ax, -2.0, 4.0)

    Colorbar(fig[1, 2i], plt, ticklabelsize=16)
    Label(fig[1, 2i, Top()], L"I_{SS}", fontsize=20)
end

## save(joinpath(pic_dir, "current_masterEq_SetII_bare2LS.png"), fig)
fig

## fig S3

fig = Figure(size=(700, 350))

configs_s2 = [
    (title="Semiclassical", dir="scan_mu_nuclear_temp"),
    (title="Lindblad", dir=joinpath(pwd(), "masterEq2", "Linblad", "scan_nuclear")),
]

for (i, c) in enumerate(configs_s2)
    ax = Axis(fig[1, 2i-1],
        xlabel = L"\mu_R/\Gamma",
        xticks = LinearTicks(4), title=c.title)
    if i == 1
        ax.ylabel = L"\lambda /\Gamma"
    end

    kwargs = isnothing(c.dir) ? () : (target_dir=c.dir,)
    _, plt = plot_steady_state_only!(ax; set="SetIII", kwargs...)

    if i == 1
        text!(ax, 0.2, 0.1, text="Set III", space=:relative, color="white")
    end
    xlims!(ax, -2.0, 4.0)

    Colorbar(fig[1, 2i], plt, ticklabelsize=16)
    Label(fig[1, 2i, Top()], L"I_{SS}", fontsize=20)
end

# save(joinpath(pic_dir, "current_masterEq_SetIII_withNuc.png"), fig)
fig

## Fig S4

function compareGauge(savefig=false)
    data_dir = "/Users/kritanjanpolley/scratch/photonFerm/withPhoton"

    f1 = loadtxt(joinpath(data_dir, "population.HEOM.test2.txt"))
    f1a = loadtxt(joinpath(data_dir, "population.HEOM.PZW_test.txt"))
    f2 = loadtxt(joinpath(data_dir, "photon.HEOM.test2.txt"))
    f2a = loadtxt(joinpath(data_dir, "photon.HEOM.PZW_test.txt"))
    f3 = loadtxt(joinpath(data_dir, "current.HEOM.test2.txt"))
    f3a = loadtxt(joinpath(data_dir, "current.HEOM.PZW_test.txt"))
    f4 = loadtxt(joinpath(data_dir, "boson_diag.HEOM.test2.txt"))
    f4a = loadtxt(joinpath(data_dir, "boson_diag.HEOM.PZW_test.txt"))

    fig = Figure(linewidth=2.2, size=(800, 300))

    ax1 = Axis(fig[1,1], 
        ylabel="Population",
        xlabel=L"\Gamma t",)
    lines!(ax1, f1[:,1], f1[:,2], label=L"\rho_{11}")
    lines!(ax1, f1[:,1], f1[:,3], label=L"\rho_{22}")
    lines!(ax1, f1a[:,1], f1a[:,2], linestyle=:dash)
    lines!(ax1, f1a[:,1], f1a[:,3], linestyle=:dash)
    axislegend(ax1)

    ax2 = Axis(fig[1,2], 
        ylabel="Photon Number",
        xlabel=L"\Gamma t",)
    lines!(ax2, f2[:,1], f2[:,2])
    lines!(ax2, f2a[:,1], f2a[:,2], linestyle=:dash)

    ax3 = Axis(fig[1,3], 
        ylabel="Current", 
        xlabel=L"\Gamma t",)
    lines!(ax3, f3[:,1], f3[:,2], label=L"I_L(t)")
    lines!(ax3, f3[:,1], f3[:,3], label=L"I_R(t)")
    lines!(ax3, f3a[:,1], f3a[:,2], linestyle=:dash)
    lines!(ax3, f3a[:,1], f3a[:,3], linestyle=:dash)
    axislegend(ax3)

    linkxaxes!(ax1, ax2, ax3)

    rowgap!(fig.layout, 2)
    colgap!(fig.layout, 4)

    # figname = "compareGauge.pdf"
    # figname = joinpath(pic_dir, figname)
    # save(figname, fig)

    
    fig
end

compareGauge()
