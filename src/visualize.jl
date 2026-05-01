using Plots
using Printf
using CSV
using JSON
using DataFrames
using Statistics

include("network.jl")
include("solvers.jl")

# ============================================================================
# Utility: Load RR placements from experiment results
# ============================================================================

"""
    load_rr_placements(results_dir) -> Dict or empty dict if not found
    
Load saved RR node placements from experiment results.
Returns empty dict if file doesn't exist.
"""
function load_rr_placements(results_dir::String)::Dict
    outpath = joinpath(results_dir, "rr_placements.json")
    if !isfile(outpath)
        return Dict()
    end
    return JSON.parse(read(outpath, String))
end

function load_topology(data_dir::String, filename::String; kwargs...)
    filepath = joinpath(data_dir, filename)
    isfile(filepath) || error("File not found: $filepath")
    gml = parse_gml(filepath)
    return build_network(gml; kwargs...)
end

# ============================================================================
# GEOGRAPHIC TOPOLOGY PLOTS
# ============================================================================

function plot_topology(net::ISPNetwork, rr_nodes::Vector{Int};
    title::String="", save_path=nothing)
    isempty(title) && (title = "$(net.name) — $(length(rr_nodes)) Route Reflectors")

    assignments = get_assignments(net, rr_nodes)
    n_rrs = length(rr_nodes)
    base_colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta,
        :brown, :pink, :olive, :teal, :navy]
    cluster_colors = base_colors[1:min(n_rrs, length(base_colors))]
    while length(cluster_colors) < n_rrs
        push!(cluster_colors, :gray)
    end
    rr_idx_map = Dict(rr => i for (i, rr) in enumerate(rr_nodes))

    p = plot(size=(1000, 700), title=title, xlabel="Longitude", ylabel="Latitude",
        legend=:outertopright, margin=5Plots.mm, dpi=200, bg=:white)

    for e in edges(net.graph)
        u, v = src(e), dst(e)
        plot!(p, [net.longitudes[u], net.longitudes[v]],
            [net.latitudes[u], net.latitudes[v]],
            color=:gray85, linewidth=0.7, label="", alpha=0.5)
    end

    for ci in 1:n_rrs
        rr = rr_nodes[ci]
        clients = [i for i in 1:net.n if assignments[i] == rr && i != rr]
        isempty(clients) && continue
        scatter!(p, net.longitudes[clients], net.latitudes[clients],
            color=cluster_colors[ci], marker=:circle, markersize=3,
            markerstrokewidth=0.2, alpha=0.6, label="")
    end

    for i in net.exit_nodes
        i in rr_nodes && continue
        ci = get(rr_idx_map, assignments[i], 1)
        scatter!(p, [net.longitudes[i]], [net.latitudes[i]],
            color=cluster_colors[ci], marker=:diamond, markersize=5,
            markerstrokewidth=0.5, alpha=0.8, label="")
    end

    for (ci, rr) in enumerate(rr_nodes)
        scatter!(p, [net.longitudes[rr]], [net.latitudes[rr]],
            color=cluster_colors[ci], marker=:star5, markersize=14,
            markerstrokecolor=:black, markerstrokewidth=2,
            label="RR: $(net.labels[rr])")
    end

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 1: Optimality Gap — histogram (capped)
# ============================================================================

function plot_gap_histogram(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment1_optimality_gap.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    greedy_gaps = min.(filter(r -> r.method == "Greedy" && !isnan(r.gap_pct), df).gap_pct, 200.0)
    ls_gaps = min.(filter(r -> r.method == "LocalSearch" && !isnan(r.gap_pct), df).gap_pct, 200.0)

    p = plot(size=(900, 550), dpi=200, margin=8Plots.mm,
        title="How Far Are Heuristics From Optimal? ($(length(greedy_gaps)) test cases)")
    histogram!(p, greedy_gaps, bins=40, alpha=0.65, label="Greedy", color=:orange)
    histogram!(p, ls_gaps, bins=40, alpha=0.65, label="Local Search", color=:green)
    xlabel!(p, "% Gap From ILP Optimal (capped at 200%)")
    ylabel!(p, "Number of Test Cases")
    vline!(p, [0], color=:black, linewidth=2, linestyle=:dash, label="Optimal (0%)")

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 1: ILP time vs size — why heuristics are needed
# ============================================================================

function plot_ilp_time(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment1_optimality_gap.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    ilp_df = filter(r -> r.method == "ILP" && !isnan(r.time_sec) && r.status != "Skipped", df)
    isempty(ilp_df) && return nothing

    p = plot(size=(900, 550), dpi=200, margin=8Plots.mm,
        title="ILP Solver Time vs Network Size (Why Heuristics Are Needed)")
    scatter!(p, ilp_df.nodes, ilp_df.time_sec,
        color=:blue, marker=:circle, markersize=4, alpha=0.5, label="ILP solve time")
    xlabel!(p, "Network Size (nodes)")
    ylabel!(p, "ILP Solve Time (seconds)")

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 2: Runtime scaling — all methods
# ============================================================================

function plot_runtime(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment2_scaling.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    p = plot(size=(900, 550), dpi=200, margin=8Plots.mm,
        title="Runtime Scaling: Which Method Is Fastest?")
    for (method, color, marker) in [("ILP", :blue, :circle),
        ("Greedy", :orange, :square),
        ("LocalSearch", :green, :diamond)]
        sub = sort(filter(r -> r.method == method, df), :nodes)
        isempty(sub) && continue
        scatter!(p, sub.nodes, sub.time_sec, label=method,
            color=color, marker=marker, markersize=4, alpha=0.6)
    end
    xlabel!(p, "Network Size (nodes)")
    ylabel!(p, "Runtime (seconds)")

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 2: Greedy vs LocalSearch head-to-head scatter
# ============================================================================

function plot_greedy_vs_ls(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment2_scaling.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    greedy_df = filter(r -> r.method == "Greedy", df)
    ls_df = filter(r -> r.method == "LocalSearch", df)
    topos = intersect(greedy_df.topology, ls_df.topology)
    g_vals, l_vals = Float64[], Float64[]
    for t in topos
        gv = filter(r -> r.topology == t, greedy_df).objective
        lv = filter(r -> r.topology == t, ls_df).objective
        (isempty(gv) || isempty(lv)) && continue
        push!(g_vals, gv[1])
        push!(l_vals, lv[1])
    end

    p = plot(size=(700, 700), dpi=200, margin=8Plots.mm,
        title="Greedy vs Local Search (each dot = 1 topology)", aspect_ratio=:equal)
    max_val = max(maximum(g_vals), maximum(l_vals)) * 1.05
    plot!(p, [0, max_val], [0, max_val], color=:red, linewidth=2, linestyle=:dash, label="Equal (y=x)")
    scatter!(p, l_vals, g_vals, color=:steelblue, markersize=4, alpha=0.5,
        label="$(length(topos)) topologies")
    xlabel!(p, "Local Search Objective (lower = better)")
    ylabel!(p, "Greedy Objective")
    annotate!(p, max_val * 0.3, max_val * 0.8, text("Points above line:\nGreedy is worse", 10, :left, :gray40))

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 3: Normalized diminishing returns
# ============================================================================

function plot_diminishing_returns(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment3_impact_of_k.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)
    topos = unique(df.topology)

    p = plot(size=(900, 550), dpi=200, margin=8Plots.mm,
        title="How Much Does Each Additional RR Help?")
    topo_colors = [:blue, :red, :green, :orange, :purple, :brown]

    for (ti, topo) in enumerate(topos)
        tsub = sort(filter(r -> r.topology == topo && r.method == "LocalSearch", df), :k)
        isempty(tsub) && continue
        base = tsub.objective[1]
        base <= 0 && continue
        pct = 100.0 .* tsub.objective ./ base
        nc = hasproperty(tsub, :nodes) && !isempty(tsub) ? tsub.nodes[1] : ""
        plot!(p, tsub.k, pct, label="$topo (n=$nc)",
            color=topo_colors[mod1(ti, length(topo_colors))],
            marker=:circle, markersize=5, linewidth=2)
    end
    xlabel!(p, "k (Number of Route Reflectors)")
    ylabel!(p, "% of k=1 Suboptimality Remaining")
    hline!(p, [0], color=:black, linewidth=1, linestyle=:dot, label="")

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 4: Dataset overview — what we tested
# ============================================================================

function plot_dataset_overview(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment4_network_structures.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    p = plot(layout=(1, 2), size=(1200, 500), dpi=200, margin=8Plots.mm)
    histogram!(p[1], df.nodes, bins=30, color=:steelblue, alpha=0.7, label="")
    xlabel!(p[1], "Network Size (nodes)")
    ylabel!(p[1], "Number of Topologies")
    title!(p[1], "Dataset: $(nrow(df)) Real ISP Topologies")
    histogram!(p[2], df.avg_degree, bins=20, color=:coral, alpha=0.7, label="")
    xlabel!(p[2], "Average Node Degree")
    ylabel!(p[2], "Number of Topologies")
    title!(p[2], "Connectivity Distribution")

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 5: Method wins — who finds best solution most often
# ============================================================================

function plot_method_wins(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment5_naive_comparison.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    methods = ["Greedy", "LocalSearch", "Betweenness", "HighDegree", "Random"]
    mc = Dict("Greedy" => :orange, "LocalSearch" => :green,
        "Betweenness" => :purple, "HighDegree" => :red, "Random" => :gray)
    k_vals = sort(unique(df.k))

    p = plot(layout=(1, length(k_vals)), size=(450 * length(k_vals), 550), dpi=200, margin=8Plots.mm)
    for (ki, k) in enumerate(k_vals)
        ksub = filter(r -> r.k == k, df)
        topos = unique(ksub.topology)
        win_counts = Dict(m => 0 for m in methods)
        for t in topos
            tsub = filter(r -> r.topology == t, ksub)
            isempty(tsub) && continue
            best_obj = minimum(tsub.objective)
            for row in eachrow(tsub)
                if row.method in methods && abs(row.objective - best_obj) < 1e-6
                    win_counts[row.method] += 1
                end
            end
        end
        counts = [get(win_counts, m, 0) for m in methods]
        colors_out = [get(mc, m, :gray) for m in methods]
        bar!(p[ki], methods, counts, color=colors_out, label="", xrotation=25)
        title!(p[ki], "k=$k: Who Finds Best? ($(length(topos)) topologies)")
        ylabel!(p[ki], "# Times Best Solution Found")
    end

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 5: Distribution strip plot — full quality spread per method
# ============================================================================

function plot_distribution(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment5_naive_comparison.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    sub = filter(r -> r.k == 3, df)
    isempty(sub) && (sub = df)
    methods = ["Greedy", "LocalSearch", "Betweenness", "HighDegree", "Random"]
    mc = Dict("Greedy" => :orange, "LocalSearch" => :green,
        "Betweenness" => :purple, "HighDegree" => :red, "Random" => :gray)

    p = plot(size=(1000, 600), dpi=200, margin=8Plots.mm,
        title="Method Quality Distribution (k=3, each dot = 1 topology)")
    for (mi, method) in enumerate(methods)
        msub = filter(r -> r.method == method, sub)
        isempty(msub) && continue
        vals = min.(msub.pct_above_best, 500.0)
        x_pos = fill(mi, length(vals)) .+ (rand(length(vals)) .- 0.5) .* 0.3
        scatter!(p, x_pos, vals, color=get(mc, method, :gray),
            markersize=3, alpha=0.4, label="")
        med = median(vals)
        plot!(p, [mi - 0.3, mi + 0.3], [med, med], color=:black, linewidth=3,
            label=(mi == 1 ? "Median" : ""))
    end
    xticks!(p, 1:length(methods), methods)
    xlabel!(p, "Method")
    ylabel!(p, "% Above Best Solution (capped at 500%)")

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# EXP 5: Before vs After — HighDegree (naive) vs Optimized
# ============================================================================

function plot_before_after(results_dir; save_path=nothing)
    csv_path = joinpath(results_dir, "experiment5_naive_comparison.csv")
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)

    sub3 = filter(r -> r.k == 3, df)
    isempty(sub3) && return nothing
    topos = unique(sub3.topology)

    naive_vals, opt_vals = Float64[], Float64[]
    for t in topos
        tsub = filter(r -> r.topology == t, sub3)
        hd = filter(r -> r.method == "HighDegree", tsub)
        ls = filter(r -> r.method == "LocalSearch", tsub)
        (isempty(hd) || isempty(ls)) && continue
        push!(naive_vals, hd.objective[1])
        push!(opt_vals, ls.objective[1])
    end

    p = plot(layout=(1, 2), size=(1400, 600), dpi=200, margin=8Plots.mm)

    if !isempty(naive_vals)
        max_val = max(maximum(naive_vals), maximum(opt_vals)) * 1.05
        plot!(p[1], [0, max_val], [0, max_val], color=:gray, linewidth=2,
            linestyle=:dash, label="No improvement")
        scatter!(p[1], opt_vals, naive_vals, color=:red, markersize=4, alpha=0.4,
            label="$(length(topos)) topologies")
        xlabel!(p[1], "Optimized (Local Search)")
        ylabel!(p[1], "Naive (Highest Degree)")
        title!(p[1], "Before vs After Optimization (k=3)")
        annotate!(p[1], max_val * 0.2, max_val * 0.85,
            text("Above line = optimization helped", 10, :left, :gray40))
    end

    improvements = Float64[]
    for (n, o) in zip(naive_vals, opt_vals)
        if o > 1e-6
            push!(improvements, 100.0 * (n - o) / n)
        elseif n > 1e-6
            push!(improvements, 100.0)
        end
    end

    if !isempty(improvements)
        histogram!(p[2], improvements, bins=30, color=:steelblue, alpha=0.7, label="")
        xlabel!(p[2], "% Improvement Over Naive HighDegree")
        ylabel!(p[2], "Number of Topologies")
        title!(p[2], "How Much Does Optimization Help?")
        vline!(p[2], [median(improvements)], color=:red, linewidth=3, linestyle=:dash,
            label="Median: $(round(median(improvements), digits=1))%")
    end

    save_path !== nothing && (savefig(p, save_path); println("  Saved: $save_path"))
    return p
end

# ============================================================================
# MASTER — Generate all 12 plots
# ============================================================================

function generate_all_plots(data_dir::String, results_dir::String, plots_dir::String)
    mkpath(plots_dir)

    println("\n" * "="^70)
    println("GENERATING VISUALIZATIONS (12 plots)")
    println("="^70)

    # Load pre-computed RR placements from experiments
    rr_placements = load_rr_placements(results_dir)

    # Check that experiments have been run
    if isempty(rr_placements)
        error("""
        ERROR: No RR placements found!

        Please run experiments FIRST to generate results:
            julia run_experiments.jl

        Then visualize the results:
            julia run_visualizations.jl

        This ensures plots are generated from actual computed data, not dummy values.
        """)
    end

    # 3 Geographic maps — using ONLY saved experimental results
    println("\n--- Geographic Maps (from Experiments) ---")
    geo_topos = [("Bellcanada.gml", 3), ("Canerie.gml", 3), ("Cogentco.gml", 5)]
    for (topo_file, k) in geo_topos
        filepath = joinpath(data_dir, topo_file)
        isfile(filepath) || continue

        name = replace(topo_file, ".gml" => "")
        net = load_topology(data_dir, topo_file)
        actual_k = min(k, net.n - 1)

        # Load RR placement from experiments (required)
        key = "$(name)_k$(actual_k)"
        if !haskey(rr_placements, key)
            error("Missing RR placement for $topo_file with k=$actual_k. Re-run experiments.")
        end

        placement = rr_placements[key]
        rr_nodes = Int.(placement["rr_nodes"])
        method_used = get(placement, "method", "LocalSearch")
        objective = placement["objective"]

        println("  Plotting $topo_file (k=$actual_k) — $(method_used) [obj=$(round(objective, digits=2))]...")
        plot_topology(net, rr_nodes;
            title="$(net.name) — $(method_used) RR Placement (k=$actual_k) [obj=$(round(objective, digits=1))]",
            save_path=joinpath(plots_dir, "geo_$(name)_k$(actual_k).png"))
    end

    # 9 Experiment plots
    println("\n--- Experiment Plots ---")
    plot_gap_histogram(results_dir; save_path=joinpath(plots_dir, "exp1_gap_histogram.png"))
    plot_ilp_time(results_dir; save_path=joinpath(plots_dir, "exp1_ilp_time.png"))
    plot_runtime(results_dir; save_path=joinpath(plots_dir, "exp2_runtime.png"))
    plot_greedy_vs_ls(results_dir; save_path=joinpath(plots_dir, "exp2_greedy_vs_ls.png"))
    plot_diminishing_returns(results_dir; save_path=joinpath(plots_dir, "exp3_diminishing_returns.png"))
    plot_dataset_overview(results_dir; save_path=joinpath(plots_dir, "exp4_dataset_overview.png"))
    plot_method_wins(results_dir; save_path=joinpath(plots_dir, "exp5_method_wins.png"))
    plot_distribution(results_dir; save_path=joinpath(plots_dir, "exp5_distribution.png"))
    plot_before_after(results_dir; save_path=joinpath(plots_dir, "exp5_before_after.png"))

    # Summary
    println("\n--- Summary ---")
    println("\n" * "="^70)
    println("DONE: 12 plots saved to $plots_dir")
    println("="^70)
end