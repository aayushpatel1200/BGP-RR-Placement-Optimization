using Pkg;
Pkg.activate(@__DIR__);
for p in ["JuMP", "HiGHS", "Graphs", "Plots", "CSV", "DataFrames", "JSON"]
    haskey(Pkg.project().dependencies, p) || Pkg.add(p)
end
# Parse command-line arguments
data_dir = length(ARGS) >= 1 ? ARGS[1] : "data"
results_dir = length(ARGS) >= 2 ? ARGS[2] : "results"

println("="^70)
println("BGP Route Reflector Optimal Placement — Experiment Runner")
println("="^70)
println("Data directory:    $data_dir")
println("Results directory: $results_dir")
println()

# Check data directory exists and has GML files
if !isdir(data_dir)
    error("""
    Data directory '$data_dir' not found!

    Please create it and place your Topology Zoo .gml files there.
    For example:
        mkdir data
        # Copy/extract your .gml files into data/

    Required topologies (minimum):
        Abilene.gml, Geant2012.gml, Dfn.gml, Deltacom.gml, Cogentco.gml, Kdl.gml

    These are available from: https://topology-zoo.org/dataset.html
    """)
end

gml_files = filter(f -> endswith(lowercase(f), ".gml"), readdir(data_dir))
if isempty(gml_files)
    error("No .gml files found in '$data_dir'. Please place Topology Zoo GML files there.")
end
println("Found $(length(gml_files)) GML files in $data_dir")

# Run experiments
cd(@__DIR__)
include("src/experiments.jl")

@time run_all_experiments(data_dir, results_dir)

println("\n\nDone! Results are in: $results_dir/")
println("Next step: run 'julia run_visualizations.jl' to generate plots.")
