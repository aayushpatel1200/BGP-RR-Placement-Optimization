"""
    solvers.jl — Three solution methods for the p-median BGP RR placement problem.
    
    Method 1: Integer Linear Program (ILP) using JuMP + HiGHS — exact solution
    Method 2: Greedy Heuristic — fast, forward selection
    Method 3: Local Search (Swap Heuristic) — iterative improvement
"""

using JuMP
using HiGHS
using Printf

# ============================================================================
# Method 1: Integer Linear Program (ILP) — Exact Solution
# ============================================================================

"""
    SolverResult — Standardized output from any solver.
"""
struct SolverResult
    method::String
    rr_nodes::Vector{Int}
    objective::Float64
    time_seconds::Float64
    status::String
end

function solve_ilp(net::ISPNetwork, k::Int; time_limit::Float64=300.0)
    n = net.n
    C = copy(net.cost_matrix)
    replace!(C, NaN => 0.0, Inf => 0.0, -Inf => 0.0)

    t_start = time()

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    set_time_limit_sec(model, time_limit)

    # Decision variables
    @variable(model, y[1:n], Bin)
    @variable(model, 0 <= x[1:n, 1:n] <= 1)
    @objective(model, Min, sum(C[i, j] * x[i, j] for i in 1:n, j in 1:n))
    @constraint(model, sum(y[j] for j in 1:n) == k)
    @constraint(model, [i = 1:n], sum(x[i, j] for j in 1:n) == 1)
    @constraint(model, [i = 1:n, j = 1:n], x[i, j] <= y[j])
    optimize!(model)

    t_elapsed = time() - t_start

    # Extract solution
    status_code = termination_status(model)
    if status_code == MOI.OPTIMAL || status_code == MOI.LOCALLY_SOLVED
        status = "Optimal"
    elseif status_code == MOI.TIME_LIMIT && has_values(model)
        status = "Feasible (timeout)"
    elseif has_values(model)
        status = "Feasible"
    else
        status = "Infeasible/Error"
        return SolverResult("ILP", Int[], Inf, t_elapsed, status)
    end

    # Get selected RR nodes
    rr_nodes = [j for j in 1:n if value(y[j]) > 0.5]
    obj = objective_value(model)

    return SolverResult("ILP", rr_nodes, obj, t_elapsed, status)
end


# ============================================================================
# Method 2: Greedy Heuristic — Forward Selection
# ============================================================================
function solve_greedy(net::ISPNetwork, k::Int)
    n = net.n
    C = copy(net.cost_matrix)
    replace!(C, NaN => 0.0, Inf => 0.0, -Inf => 0.0)

    t_start = time()

    selected = Int[]
    remaining = Set(1:n)
    client_cost = fill(Inf, n)

    for step in 1:k
        best_node = 0

        if step == 1
            best_total = Inf
            for j in remaining
                total = 0.0
                for i in 1:n
                    total += C[i, j]
                end
                if total < best_total
                    best_total = total
                    best_node = j
                end
            end
        else
            best_reduction = -Inf
            for j in remaining
                reduction = 0.0
                for i in 1:n
                    new_cost = min(client_cost[i], C[i, j])
                    reduction += (client_cost[i] - new_cost)
                end
                if reduction > best_reduction
                    best_reduction = reduction
                    best_node = j
                end
            end
        end
        push!(selected, best_node)
        delete!(remaining, best_node)
        for i in 1:n
            client_cost[i] = min(client_cost[i], C[i, best_node])
        end
    end

    t_elapsed = time() - t_start
    obj = sum(client_cost)

    return SolverResult("Greedy", selected, obj, t_elapsed, "Heuristic")
end


# ============================================================================
# Method 3: Local Search (Swap Heuristic)
# ============================================================================
function solve_local_search(net::ISPNetwork, k::Int;
    max_iter::Int=1000, restarts::Int=5, seed::Int=42)
    rng = MersenneTwister(seed)
    n = net.n
    C = copy(net.cost_matrix)
    replace!(C, NaN => 0.0, Inf => 0.0, -Inf => 0.0)

    t_start = time()

    best_global_rrs = Int[]
    best_global_obj = Inf

    for restart in 1:restarts
        current_rrs = sort(shuffle(rng, collect(1:n))[1:k])
        current_obj = evaluate_placement(net, current_rrs)

        improved = true
        iter = 0

        while improved && iter < max_iter
            improved = false
            iter += 1

            best_swap_gain = 0.0
            best_add = 0
            best_remove = 0

            rr_set = Set(current_rrs)
            non_rrs = [v for v in 1:n if !(v in rr_set)]

            for (ri, r) in enumerate(current_rrs)
                for a in non_rrs
                    new_rrs = copy(current_rrs)
                    new_rrs[ri] = a
                    new_obj = evaluate_placement(net, new_rrs)

                    gain = current_obj - new_obj
                    if gain > best_swap_gain
                        best_swap_gain = gain
                        best_remove = ri
                        best_add = a
                    end
                end
            end

            if best_swap_gain > 1e-10
                current_rrs[best_remove] = best_add
                sort!(current_rrs)
                current_obj -= best_swap_gain
                current_obj = evaluate_placement(net, current_rrs)
                improved = true
            end
        end

        if current_obj < best_global_obj
            best_global_obj = current_obj
            best_global_rrs = copy(current_rrs)
        end
    end

    t_elapsed = time() - t_start
    final_obj = evaluate_placement(net, best_global_rrs)

    return SolverResult("LocalSearch", best_global_rrs, final_obj, t_elapsed, "Heuristic")
end


# ============================================================================
# Naive Baselines (for Experiment 5)
# ============================================================================
function solve_highest_betweenness(net::ISPNetwork, k::Int)
    t_start = time()

    bc = betweenness_centrality(net.graph)
    sorted_nodes = sortperm(bc, rev=true)
    rr_nodes = sorted_nodes[1:min(k, length(sorted_nodes))]

    obj = evaluate_placement(net, rr_nodes)
    t_elapsed = time() - t_start

    return SolverResult("Betweenness", rr_nodes, obj, t_elapsed, "Baseline")
end

function solve_highest_degree(net::ISPNetwork, k::Int)
    t_start = time()

    deg = [degree(net.graph, v) for v in 1:net.n]
    sorted_nodes = sortperm(deg, rev=true)
    rr_nodes = sorted_nodes[1:min(k, length(sorted_nodes))]

    obj = evaluate_placement(net, rr_nodes)
    t_elapsed = time() - t_start

    return SolverResult("HighDegree", rr_nodes, obj, t_elapsed, "Baseline")
end

function solve_random(net::ISPNetwork, k::Int; num_trials::Int=20, seed::Int=123)
    rng = MersenneTwister(seed)
    t_start = time()

    best_obj = Inf
    best_rrs = Int[]

    for _ in 1:num_trials
        rrs = sort(shuffle(rng, collect(1:net.n))[1:k])
        obj = evaluate_placement(net, rrs)
        if obj < best_obj
            best_obj = obj
            best_rrs = rrs
        end
    end

    t_elapsed = time() - t_start
    return SolverResult("Random", best_rrs, best_obj, t_elapsed, "Baseline")
end
