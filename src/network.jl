"""
    network.jl — Network model for BGP Route Reflector placement.
    
    Builds a weighted graph from GML data, computes all-pairs shortest paths,
    identifies exit points, and computes the routing suboptimality cost matrix.
"""

using Graphs
using Printf
using Random
using Statistics

include("gml_parser.jl")

struct ISPNetwork
    name::String
    n::Int
    labels::Vector{String}
    latitudes::Vector{Float64}
    longitudes::Vector{Float64}
    graph::SimpleGraph{Int}
    weights::Matrix{Float64}
    dist::Matrix{Float64}
    exit_nodes::Vector{Int}
    cost_matrix::Matrix{Float64}
end

function haversine(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)
    R = 6371.0  # Earth radius in km
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    a = sin(dlat / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return R * c
end

function build_network(gml::GMLGraph; exit_fraction::Float64=0.20, min_exits::Int=3,
    max_prefix_pairs::Int=500, seed::Int=42)
    rng = MersenneTwister(seed)
    n = length(gml.nodes)
    id_to_idx = Dict{Int,Int}()
    labels = String[]
    lats = Float64[]
    lons = Float64[]

    for (idx, node) in enumerate(gml.nodes)
        id_to_idx[node.id] = idx
        push!(labels, node.label)
        push!(lats, node.latitude)
        push!(lons, node.longitude)
    end
    valid_lats = filter(!isnan, lats)
    valid_lons = filter(!isnan, lons)
    mean_lat = isempty(valid_lats) ? 0.0 : mean(valid_lats)
    mean_lon = isempty(valid_lons) ? 0.0 : mean(valid_lons)
    for i in 1:n
        if isnan(lats[i])
            lats[i] = mean_lat + randn(rng) * 0.1
            lons[i] = mean_lon + randn(rng) * 0.1
        end
    end

    g = SimpleGraph(n)
    weight_matrix = fill(Inf, n, n)
    for i in 1:n
        weight_matrix[i, i] = 0.0
    end

    for edge in gml.edges
        if !haskey(id_to_idx, edge.source) || !haskey(id_to_idx, edge.target)
            continue
        end
        u = id_to_idx[edge.source]
        v = id_to_idx[edge.target]
        if u == v
            continue
        end
        add_edge!(g, u, v)
        w = haversine(lats[u], lons[u], lats[v], lons[v])
        w = max(w, 1.0)
        weight_matrix[u, v] = w
        weight_matrix[v, u] = w
    end

    dist_matrix = fill(Inf, n, n)
    for s in 1:n
        dists = dijkstra_custom(n, s, weight_matrix, g)
        dist_matrix[s, :] = dists
    end

    deg1_nodes = [v for v in 1:n if degree(g, v) == 1]

    target_exits = max(min_exits, ceil(Int, exit_fraction * n))

    if length(deg1_nodes) >= min_exits
        exit_nodes = deg1_nodes
    else
        candidates = setdiff(1:n, deg1_nodes)
        num_extra = target_exits - length(deg1_nodes)
        num_extra = min(num_extra, length(candidates))
        extra = shuffle(rng, collect(candidates))[1:num_extra]
        exit_nodes = sort(vcat(deg1_nodes, extra))
    end

    if length(exit_nodes) > 50
        exit_nodes = sort(shuffle(rng, exit_nodes)[1:50])
    end

    @printf("  Network: %s | Nodes: %d | Edges: %d | Exits: %d\n",
        gml.name, n, ne(g), length(exit_nodes))

    cost_mat = compute_cost_matrix(n, dist_matrix, exit_nodes, max_prefix_pairs, rng)

    return ISPNetwork(gml.name, n, labels, lats, lons, g, weight_matrix,
        dist_matrix, exit_nodes, cost_mat)
end

function dijkstra_custom(n::Int, src::Int, weight_matrix::Matrix{Float64}, g::SimpleGraph)
    dist = fill(Inf, n)
    visited = falses(n)
    dist[src] = 0.0

    for _ in 1:n
        u = 0
        min_d = Inf
        for v in 1:n
            if !visited[v] && dist[v] < min_d
                min_d = dist[v]
                u = v
            end
        end
        if u == 0 || min_d == Inf
            break
        end
        visited[u] = true

        for v in neighbors(g, u)
            if !visited[v]
                new_dist = dist[u] + weight_matrix[u, v]
                if new_dist < dist[v]
                    dist[v] = new_dist
                end
            end
        end
    end

    return dist
end

function compute_cost_matrix(n::Int, dist::Matrix{Float64}, exits::Vector{Int},
    max_pairs::Int, rng::AbstractRNG)
    cost = zeros(Float64, n, n)
    ne = length(exits)

    if ne < 2
        # Fallback: simple distance-based cost
        if ne == 1
            for i in 1:n, j in 1:n
                cost[i, j] = 0.0
            end
        end
        return cost
    end

    all_pairs = [(exits[a], exits[b]) for a in 1:ne for b in (a+1):ne]

    if length(all_pairs) > max_pairs
        all_pairs = shuffle(rng, all_pairs)[1:max_pairs]
    end

    num_pairs = length(all_pairs)

    for (ea, eb) in all_pairs
        for j in 1:n  # RR location
            if isinf(dist[j, ea]) && isinf(dist[j, eb])
                continue
            end

            rr_picks_ea = dist[j, ea] <= dist[j, eb]

            for i in 1:n  # Client location
                d_i_ea = dist[i, ea]
                d_i_eb = dist[i, eb]
                if isinf(d_i_ea) && isinf(d_i_eb)
                    continue
                end

                client_optimal = min(d_i_ea, d_i_eb)
                rr_exit_cost = rr_picks_ea ? d_i_ea : d_i_eb
                if isinf(rr_exit_cost)
                    continue
                end
                penalty = max(0.0, rr_exit_cost - client_optimal)
                cost[i, j] += penalty
            end
        end
    end

    if num_pairs > 0
        cost ./= num_pairs
    end

    return cost
end

function evaluate_placement(net::ISPNetwork, rr_nodes::Vector{Int})
    total_cost = 0.0
    for i in 1:net.n
        best_cost = Inf
        for j in rr_nodes
            c = net.cost_matrix[i, j]
            if c < best_cost
                best_cost = c
            end
        end
        total_cost += best_cost
    end
    return total_cost
end

function get_assignments(net::ISPNetwork, rr_nodes::Vector{Int})
    assignments = zeros(Int, net.n)
    for i in 1:net.n
        best_cost = Inf
        best_rr = rr_nodes[1]
        for j in rr_nodes
            c = net.cost_matrix[i, j]
            if c < best_cost
                best_cost = c
                best_rr = j
            end
        end
        assignments[i] = best_rr
    end
    return assignments
end
