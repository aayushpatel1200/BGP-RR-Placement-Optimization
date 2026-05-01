"""
    experiments.jl — Run all 5 experiments for the BGP RR placement project.
    
    Experiment 1: Optimality Gap (ILP vs Greedy vs Local Search)
    Experiment 2: Scaling Behavior  
    Experiment 3: Impact of k (number of RRs)
    Experiment 4: Different Network Structures
    Experiment 5: Comparison to Naive Placement
"""

using Printf
using CSV
using JSON
using DataFrames
using Statistics

include("network.jl")
include("solvers.jl")

# ============================================================================
# Utility: Save RR placements for visualization
# ============================================================================
function save_rr_placements(results_dir::String, placements::Dict)
    mkpath(results_dir)
    outpath = joinpath(results_dir, "rr_placements.json")
    open(outpath, "w") do f
        write(f, JSON.json(placements, 2))
    end
end

function load_rr_placements(results_dir::String)::Dict
    outpath = joinpath(results_dir, "rr_placements.json")
    if !isfile(outpath)
        return Dict()
    end
    return JSON.parse(read(outpath, String))
end

# ============================================================================
# Utility: Load a topology by filename
# ============================================================================
function load_topology(data_dir::String, filename::String; kwargs...)
    filepath = joinpath(data_dir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    gml = parse_gml(filepath)
    return build_network(gml; kwargs...)
end

# ============================================================================
# Experiment 1: Optimality Gap — ILP vs Greedy vs Local Search
# ============================================================================
function experiment1_optimality_gap(data_dir::String, results_dir::String;
    topologies=["Abilene.gml", "Geant2012.gml", "Dfn.gml", "Deltacom.gml"],
    k_values=[1, 2, 3, 5],
    ilp_time_limit=120.0)

    println("\n" * "="^70)
    println("EXPERIMENT 1: Optimality Gap — ILP vs Greedy vs Local Search")
    println("="^70)

    rows = []

    for topo_file in topologies
        filepath = joinpath(data_dir, topo_file)
        if !isfile(filepath)
            println("  WARNING: Skipping $topo_file (file not found)")
            continue
        end

        println("\nLoading $topo_file...")
        net = load_topology(data_dir, topo_file)

        for k in k_values
            if k >= net.n
                continue
            end

            println("  k=$k:")

            ilp_result = nothing
            if net.n <= 150 && k <= 10
                print("    ILP... ")
                ilp_result = solve_ilp(net, k; time_limit=ilp_time_limit)
                @printf("obj=%.2f  time=%.2fs  status=%s\n",
                    ilp_result.objective, ilp_result.time_seconds, ilp_result.status)
            else
                println("    ILP... SKIPPED (n=$(net.n) too large)")
            end

            # Greedy
            print("    Greedy... ")
            greedy_result = solve_greedy(net, k)
            @printf("obj=%.2f  time=%.2fs\n", greedy_result.objective, greedy_result.time_seconds)

            # Local Search
            print("    Local Search... ")
            ls_restarts = net.n <= 200 ? 5 : 3
            ls_result = solve_local_search(net, k; restarts=ls_restarts)
            @printf("obj=%.2f  time=%.2fs\n", ls_result.objective, ls_result.time_seconds)

            # Compute gap relative to best known
            best_obj = min(greedy_result.objective, ls_result.objective)
            if ilp_result !== nothing && ilp_result.objective < Inf
                best_obj = min(best_obj, ilp_result.objective)
            end

            # Record results
            for result in [ilp_result, greedy_result, ls_result]
                if result === nothing
                    push!(rows, (topology=net.name, nodes=net.n, k=k,
                        method="ILP", objective=NaN, time_sec=NaN,
                        gap_pct=NaN, status="Skipped"))
                    continue
                end
                gap = best_obj > 1e-10 ?
                      100.0 * (result.objective - best_obj) / best_obj : 0.0
                push!(rows, (topology=net.name, nodes=net.n, k=k,
                    method=result.method, objective=result.objective,
                    time_sec=result.time_seconds, gap_pct=gap,
                    status=result.status))
            end
        end
    end

    df = DataFrame(rows)
    outpath = joinpath(results_dir, "experiment1_optimality_gap.csv")
    CSV.write(outpath, df)
    println("\nResults saved to $outpath")
    return df
end


# ============================================================================
# Experiment 2: Scaling Behavior
# ============================================================================
function experiment2_scaling(data_dir::String, results_dir::String;
    topologies=["Abilene.gml", "Geant2012.gml", "Dfn.gml",
        "Deltacom.gml", "Cogentco.gml", "Kdl.gml"],
    k=3)

    println("\n" * "="^70)
    println("EXPERIMENT 2: Scaling Behavior")
    println("="^70)

    rows = []

    for topo_file in topologies
        filepath = joinpath(data_dir, topo_file)
        if !isfile(filepath)
            println("  WARNING: Skipping $topo_file (file not found)")
            continue
        end

        println("\nLoading $topo_file...")
        net = load_topology(data_dir, topo_file)

        actual_k = min(k, net.n - 1)

        # Greedy
        print("  Greedy (k=$actual_k)... ")
        greedy_result = solve_greedy(net, actual_k)
        @printf("obj=%.2f  time=%.4fs\n", greedy_result.objective, greedy_result.time_seconds)
        push!(rows, (topology=net.name, nodes=net.n, edges=ne(net.graph),
            method="Greedy", k=actual_k, objective=greedy_result.objective,
            time_sec=greedy_result.time_seconds))

        # Local Search
        print("  Local Search (k=$actual_k)... ")
        ls_restarts = net.n <= 200 ? 5 : 2
        ls_result = solve_local_search(net, actual_k; restarts=ls_restarts)
        @printf("obj=%.2f  time=%.4fs\n", ls_result.objective, ls_result.time_seconds)
        push!(rows, (topology=net.name, nodes=net.n, edges=ne(net.graph),
            method="LocalSearch", k=actual_k, objective=ls_result.objective,
            time_sec=ls_result.time_seconds))

        # ILP only for small networks
        if net.n <= 150
            print("  ILP (k=$actual_k)... ")
            ilp_result = solve_ilp(net, actual_k; time_limit=180.0)
            @printf("obj=%.2f  time=%.4fs  status=%s\n",
                ilp_result.objective, ilp_result.time_seconds, ilp_result.status)
            push!(rows, (topology=net.name, nodes=net.n, edges=ne(net.graph),
                method="ILP", k=actual_k, objective=ilp_result.objective,
                time_sec=ilp_result.time_seconds))
        end
    end

    df = DataFrame(rows)
    outpath = joinpath(results_dir, "experiment2_scaling.csv")
    CSV.write(outpath, df)
    println("\nResults saved to $outpath")
    return df
end


# ============================================================================
# Experiment 3: Impact of k (Number of RRs)
# ============================================================================
function experiment3_impact_of_k(data_dir::String, results_dir::String;
    topologies=["Cogentco.gml"],
    k_range=[1, 2, 3, 5, 7, 10, 15, 20])

    println("\n" * "="^70)
    println("EXPERIMENT 3: Impact of k (Number of Route Reflectors)")
    println("="^70)

    rows = []

    for topo_file in topologies
        filepath = joinpath(data_dir, topo_file)
        if !isfile(filepath)
            println("  WARNING: Skipping $topo_file (not found)")
            continue
        end

        println("\nLoading $topo_file...")
        net = load_topology(data_dir, topo_file)

        for k in k_range
            if k >= net.n
                continue
            end

            println("  k = $k:")

            print("    Greedy... ")
            g_res = solve_greedy(net, k)
            @printf("obj=%.2f  time=%.2fs\n", g_res.objective, g_res.time_seconds)
            push!(rows, (topology=net.name, nodes=net.n, k=k, method="Greedy",
                objective=g_res.objective, time_sec=g_res.time_seconds))

            print("    Local Search... ")
            ls_res = solve_local_search(net, k; restarts=3)
            @printf("obj=%.2f  time=%.2fs\n", ls_res.objective, ls_res.time_seconds)
            push!(rows, (topology=net.name, nodes=net.n, k=k, method="LocalSearch",
                objective=ls_res.objective, time_sec=ls_res.time_seconds))

            if net.n <= 300 && k <= 10
                print("    ILP... ")
                ilp_res = solve_ilp(net, k; time_limit=120.0)
                @printf("obj=%.2f  time=%.2fs  %s\n",
                    ilp_res.objective, ilp_res.time_seconds, ilp_res.status)
                push!(rows, (topology=net.name, nodes=net.n, k=k, method="ILP",
                    objective=ilp_res.objective, time_sec=ilp_res.time_seconds))
            end
        end
    end

    df = DataFrame(rows)
    outpath = joinpath(results_dir, "experiment3_impact_of_k.csv")
    CSV.write(outpath, df)
    println("\nResults saved to $outpath")
    return df
end


# ============================================================================
# Experiment 4: Different Network Structures
# ============================================================================
function experiment4_network_structures(data_dir::String, results_dir::String;
    topologies=["Abilene.gml", "Geant2012.gml", "Dfn.gml",
        "Deltacom.gml", "Cogentco.gml", "Kdl.gml"],
    k=3)

    println("\n" * "="^70)
    println("EXPERIMENT 4: Different Network Structures")
    println("="^70)

    rows = []

    for topo_file in topologies
        filepath = joinpath(data_dir, topo_file)
        if !isfile(filepath)
            println("  WARNING: Skipping $topo_file (file not found)")
            continue
        end

        println("\nLoading $topo_file...")
        net = load_topology(data_dir, topo_file)

        actual_k = min(k, net.n - 1)

        # Compute network properties
        avg_degree = 2.0 * ne(net.graph) / net.n
        density = 2.0 * ne(net.graph) / (net.n * (net.n - 1))
        bc = betweenness_centrality(net.graph)
        max_bc = maximum(bc)

        # Greedy
        g_res = solve_greedy(net, actual_k)

        # Local Search
        ls_res = solve_local_search(net, actual_k; restarts=3)

        best_obj = min(g_res.objective, ls_res.objective)

        push!(rows, (topology=net.name, nodes=net.n, edges=ne(net.graph),
            avg_degree=round(avg_degree, digits=2),
            density=round(density, digits=4),
            max_betweenness=round(max_bc, digits=4),
            num_exits=length(net.exit_nodes),
            k=actual_k,
            greedy_obj=round(g_res.objective, digits=2),
            ls_obj=round(ls_res.objective, digits=2),
            best_obj=round(best_obj, digits=2),
            greedy_time=round(g_res.time_seconds, digits=4),
            ls_time=round(ls_res.time_seconds, digits=4)))

        @printf("  %s: n=%d, e=%d, avg_deg=%.1f, best_obj=%.2f\n",
            net.name, net.n, ne(net.graph), avg_degree, best_obj)
    end

    df = DataFrame(rows)
    outpath = joinpath(results_dir, "experiment4_network_structures.csv")
    CSV.write(outpath, df)
    println("\nResults saved to $outpath")
    return df
end


# ============================================================================
# Experiment 5: Comparison to Naive Placement
# ============================================================================
function experiment5_naive_comparison(data_dir::String, results_dir::String;
    topologies=["Abilene.gml", "Geant2012.gml", "Dfn.gml",
        "Deltacom.gml", "Cogentco.gml"],
    k_values=[1, 3, 5])

    println("\n" * "="^70)
    println("EXPERIMENT 5: Comparison to Naive Placement Heuristics")
    println("="^70)

    rows = []

    for topo_file in topologies
        filepath = joinpath(data_dir, topo_file)
        if !isfile(filepath)
            println("  WARNING: Skipping $topo_file (file not found)")
            continue
        end

        println("\nLoading $topo_file...")
        net = load_topology(data_dir, topo_file)

        for k in k_values
            if k >= net.n
                continue
            end

            println("  k=$k:")

            # Optimized methods
            greedy_res = solve_greedy(net, k)
            ls_res = solve_local_search(net, k; restarts=3)
            best_optimized = min(greedy_res.objective, ls_res.objective)

            # Naive baselines
            bc_res = solve_highest_betweenness(net, k)
            deg_res = solve_highest_degree(net, k)
            rand_res = solve_random(net, k)

            for result in [greedy_res, ls_res, bc_res, deg_res, rand_res]
                improvement = best_optimized > 1e-10 ?
                              100.0 * (result.objective - best_optimized) / best_optimized : 0.0

                push!(rows, (topology=net.name, nodes=net.n, k=k,
                    method=result.method, objective=round(result.objective, digits=2),
                    pct_above_best=round(improvement, digits=2),
                    time_sec=round(result.time_seconds, digits=4)))

                @printf("    %-15s obj=%.2f  (+%.1f%%)\n",
                    result.method, result.objective, improvement)
            end
        end
    end

    df = DataFrame(rows)
    outpath = joinpath(results_dir, "experiment5_naive_comparison.csv")
    CSV.write(outpath, df)
    println("\nResults saved to $outpath")
    return df
end


# ============================================================================
# Run All Experiments
# ============================================================================
function run_all_experiments(data_dir::String, results_dir::String)
    mkpath(results_dir)

    # Use ALL available topologies
    all_gml = sort(filter(f -> endswith(lowercase(f), ".gml"), readdir(data_dir)))
    if isempty(all_gml)
        error("No .gml files found in $data_dir. Please place Topology Zoo GML files there.")
    end

    # Sort by network size for organized output
    topo_sizes = Dict{String,Int}()
    for f in all_gml
        gml = parse_gml(joinpath(data_dir, f))
        topo_sizes[f] = length(gml.nodes)
    end
    sort!(all_gml, by=f -> topo_sizes[f])

    # Split by size for ILP feasibility
    small_topos = filter(f -> topo_sizes[f] <= 150, all_gml)   # ILP feasible
    medium_topos = filter(f -> 150 < topo_sizes[f] <= 400, all_gml)
    large_topos = filter(f -> topo_sizes[f] > 400, all_gml)

    println("\n" * "#"^70)
    println("# Starting BGP Route Reflector Placement Experiments")
    println("# Data directory: $data_dir")
    println("# Results directory: $results_dir")
    println("# Total topologies: $(length(all_gml))")
    println("# Small (≤150 nodes): $(length(small_topos))")
    println("# Medium (151-400):   $(length(medium_topos))")
    println("# Large (>400):       $(length(large_topos))")
    println("#"^70)

    # Experiment 1: Optimality Gap — ILP vs heuristics (small + medium topos)
    exp1_topos = vcat(small_topos, medium_topos)
    df1 = experiment1_optimality_gap(data_dir, results_dir;
        topologies=exp1_topos, k_values=[1, 2, 3, 5])

    # Experiment 2: Scaling — ALL topologies, heuristics only for large
    df2 = experiment2_scaling(data_dir, results_dir;
        topologies=all_gml, k=3)

    # Experiment 3: Impact of k — run on 3 diverse topologies
    exp3_topos = String[]
    for target in [40, 120, 200]
        best = argmin(f -> abs(topo_sizes[f] - target), all_gml)
        push!(exp3_topos, best)
    end
    df3 = experiment3_impact_of_k(data_dir, results_dir;
        topologies=unique(exp3_topos), k_range=[1, 2, 3, 5, 7, 10, 15, 20])

    # Experiment 4: Network structures — ALL topologies
    df4 = experiment4_network_structures(data_dir, results_dir;
        topologies=all_gml, k=3)

    # Experiment 5: Naive comparison — ALL topologies
    df5 = experiment5_naive_comparison(data_dir, results_dir;
        topologies=all_gml, k_values=[1, 3, 5])

    # Save RR placements for geographic visualizations
    println("\nGenerating RR placements for geographic visualizations...")
    placements_data = Dict()
    geo_topos = [("Bellcanada.gml", 3), ("Canerie.gml", 3), ("Cogentco.gml", 5)]
    for (topo_file, k) in geo_topos
        filepath = joinpath(data_dir, topo_file)
        if isfile(filepath)
            try
                net = load_topology(data_dir, topo_file)
                actual_k = min(k, net.n - 1)
                # Use Local Search for best quality visualization results
                result = solve_local_search(net, actual_k; restarts=5)
                key = "$(replace(topo_file, ".gml" => ""))_k$(actual_k)"
                placements_data[key] = Dict(
                    "topology" => net.name,
                    "nodes" => net.n,
                    "k" => actual_k,
                    "method" => result.method,
                    "rr_nodes" => result.rr_nodes,
                    "objective" => result.objective
                )
                println("  Saved $key (method=$(result.method))")
            catch e
                println("  WARNING: Failed to process $topo_file: $e")
            end
        end
    end
    save_rr_placements(results_dir, placements_data)

    println("\n" * "="^70)
    println("ALL EXPERIMENTS COMPLETE")
    println("Results saved in: $results_dir")
    println("="^70)

    return (df1, df2, df3, df4, df5)
end
