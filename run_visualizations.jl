using Pkg;
Pkg.activate(@__DIR__);
for p in ["JuMP", "HiGHS", "Graphs", "Plots", "CSV", "DataFrames", "JSON"]
    haskey(Pkg.project().dependencies, p) || Pkg.add(p)
end

data_dir = length(ARGS) >= 1 ? ARGS[1] : "data"
results_dir = length(ARGS) >= 2 ? ARGS[2] : "results"
plots_dir = length(ARGS) >= 3 ? ARGS[3] : joinpath(results_dir, "plots")

println("="^70)
println("BGP Route Reflector Placement — Visualization Generator")
println("="^70)
println("Data directory:    $data_dir")
println("Results directory: $results_dir")
println("Plots directory:   $plots_dir")
println()

cd(@__DIR__)
include("src/visualize.jl")

@time generate_all_plots(data_dir, results_dir, plots_dir)

println("\n\nDone! Plots saved in: $plots_dir/")
