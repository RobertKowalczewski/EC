using Random, ProgressMeter, Statistics
using Base.Threads

include("helpers.jl")
include("moves.jl")

const NUM_RUNS = 20

fmt(x) = round(x; digits=5)

"""Run independent starts in parallel and capture their statistics."""
function run_parallel_runs(num_runs, distance_matrix, costs, run_solver; record_iterations::Bool=false, start_generator=random_start)
    objectives = Vector{Float64}(undef, num_runs)
    solutions = Vector{Vector{Int}}(undef, num_runs)
    starts = Vector{Vector{Int}}(undef, num_runs)
    start_objectives = Vector{Float64}(undef, num_runs)
    times = Vector{Float64}(undef, num_runs)
    iter_counts = record_iterations ? Vector{Int}(undef, num_runs) : nothing

    progress = Progress(num_runs; desc="runs", dt=0.1)
    progress_lock = SpinLock()

    @threads for run in 1:num_runs
        rng = Random.TaskLocalRNG()
        start_solution = start_generator(distance_matrix; rng=rng)
        starts[run] = copy(start_solution)
        start_objectives[run] = calculate_cycle_length(start_solution, distance_matrix, costs)

        elapsed_time = @elapsed begin
            result = run_solver(start_solution, rng)
            if record_iterations
                objective, solution, iterations = result
                iter_counts[run] = iterations
            else
                objective, solution = result
            end
            objectives[run] = Float64(objective)
            solutions[run] = solution
        end

        times[run] = elapsed_time

        lock(progress_lock) do
            ProgressMeter.next!(progress)
        end
    end

    return (; objectives, solutions, starts, start_objectives, times, iterations=iter_counts)
end

function summarize_runs(run_data)
    objectives = run_data.objectives
    start_objectives = run_data.start_objectives
    times = run_data.times

    min_idx = argmin(objectives)
    max_idx = argmax(objectives)

    return (
        avg_objective=mean(objectives),
        min_objective=objectives[min_idx],
        max_objective=objectives[max_idx],
        avg_start=mean(start_objectives),
        min_start=minimum(start_objectives),
        max_start=maximum(start_objectives),
        avg_time=mean(times),
        min_time=minimum(times),
        max_time=maximum(times),
        best_index=min_idx,
    )
end

function test_algorithms(distance_matrix, costs, data_path, data_name)
    neighbors = precompute_nearest_neighbors(distance_matrix, 15)

    println("MSLS:")
    msls_runs = run_parallel_runs(
        NUM_RUNS,
        distance_matrix,
        costs,
        (start_solution, rng) -> begin
            multiple_start_local_search(distance_matrix, costs; rng=rng)
        end,
    )
    msls_stats = summarize_runs(msls_runs)
    report_algorithm(msls_stats)
    plot_best_solution(msls_runs.objectives, msls_runs.solutions,
        data_path,
        "$(data_name)_MSLS",
        "$(data_name)_multiple_start_LS.png")

    println("solution from first node:")
    println([x - 1 for x in msls_runs.solutions[1]])

    run_time = max(msls_stats.avg_time, eps())
    println("run time: ", run_time)

    println("ILS (Optimized):")
    ils_runs = run_parallel_runs(
        NUM_RUNS,
        distance_matrix,
        costs,
        (start_solution, rng) -> begin
            iterated_local_search(distance_matrix, costs, run_time, neighbors;
                rng=rng,
                initial_solution=copy(start_solution))
        end;
        record_iterations=true,
        start_generator=(dm; rng) -> construct_greedy_solution(dm, costs; rng=rng)
    )
    ils_stats = summarize_runs(ils_runs)
    report_algorithm(ils_stats)
    println("number of basic iterations of LS: $(ils_runs.iterations)")
    plot_best_solution(ils_runs.objectives, ils_runs.solutions,
        data_path,
        "$(data_name)_ILS",
        "$(data_name)_iterated_local_search_recomputed_time.png")

    println("ILS (Memory):")
    ils_mem_runs = run_parallel_runs(
        NUM_RUNS,
        distance_matrix,
        costs,
        (start_solution, rng) -> begin
            iterated_local_search(distance_matrix, costs, run_time, neighbors;
                rng=rng,
                initial_solution=copy(start_solution),
                use_memory=true)
        end;
        record_iterations=true,
        start_generator=(dm; rng) -> construct_greedy_solution(dm, costs; rng=rng)
    )
    ils_mem_stats = summarize_runs(ils_mem_runs)
    report_algorithm(ils_mem_stats)
    println("number of basic iterations of LS: $(ils_mem_runs.iterations)")
    plot_best_solution(ils_mem_runs.objectives, ils_mem_runs.solutions,
        data_path,
        "$(data_name)_ILS_Memory",
        "$(data_name)_iterated_local_search_memory_recomputed_time.png")

    println("best_solution")
    println([ils_runs.solutions[ils_stats.best_index]])
end

function report_algorithm(stats)
    println("After Local Search | Starting Solution")
    println("Obj: $(fmt(stats.avg_objective)) ($(fmt(stats.min_objective)) - $(fmt(stats.max_objective))) | Obj: $(fmt(stats.avg_start)) ($(fmt(stats.min_start)) - $(fmt(stats.max_start)))")
    println("Time[s]: $(fmt(stats.avg_time)) ($(fmt(stats.min_time)) - $(fmt(stats.max_time)))")
end

"""Apply structured perturbations (2-opt, swaps, replacements) to diversify the route."""
function perturb(solution, distance_matrix, perturb_size; rng=Random.default_rng())
    perturbed = copy(solution)
    n = length(solution)
    total_nodes = size(distance_matrix, 1)
    steps = max(2, 2 * Int(clamp(round(perturb_size), 1, n)))

    available_nodes = Vector{Int}()
    in_solution = falses(total_nodes)
    @inbounds for node in perturbed
        in_solution[node] = true
    end
    for node in 1:total_nodes
        if !in_solution[node]
            push!(available_nodes, node)
        end
    end

    for _ in 1:steps
        r = rand(rng)
        if r < 0.4 && n >= 4
            i = rand(rng, 1:n-2)
            j = rand(rng, i+2:n)
            reverse!(perturbed, i, j)
        elseif r < 0.7 && n >= 2
            i = rand(rng, 1:n)
            j = rand(rng, 1:n)
            while j == i
                j = rand(rng, 1:n)
            end
            perturbed[i], perturbed[j] = perturbed[j], perturbed[i]
        elseif !isempty(available_nodes)
            idx = rand(rng, 1:n)
            avail_idx = rand(rng, 1:length(available_nodes))
            new_node = available_nodes[avail_idx]
            old_node = perturbed[idx]
            perturbed[idx] = new_node
            available_nodes[avail_idx] = old_node
        end
    end

    return perturbed
end


function precompute_nearest_neighbors(distance_matrix, k)
    n = size(distance_matrix, 1)
    neighbors = Matrix{Int}(undef, n, k)
    for i in 1:n
        dists = Vector{Tuple{Float64, Int}}(undef, n-1)
        idx = 1
        for j in 1:n
            if i != j
                dists[idx] = (distance_matrix[i, j], j)
                idx += 1
            end
        end
        partialsort!(dists, 1:k, by=x->x[1])
        for m in 1:k
            neighbors[i, m] = dists[m][2]
        end
    end
    return neighbors
end

function construct_greedy_solution(distance_matrix, costs; rng=Random.default_rng())
    n = size(distance_matrix, 1)
    start_node = rand(rng, 1:n)
    solution = [start_node]
    in_solution = falses(n)
    in_solution[start_node] = true
    target_size = ceil(Int, n / 2)
    
    while length(solution) < target_size
        best_cost = Inf
        best_node = -1
        best_pos = -1
        
        for node in 1:n
            if !in_solution[node]
                for pos in 1:length(solution)
                    prev = solution[pos]
                    next = solution[mod1(pos+1, length(solution))]
                    cost_increase = distance_matrix[prev, node] + distance_matrix[node, next] - distance_matrix[prev, next] + costs[node]
                    
                    if cost_increase < best_cost
                        best_cost = cost_increase
                        best_node = node
                        best_pos = pos
                    end
                end
            end
        end
        
        insert!(solution, best_pos + 1, best_node)
        in_solution[best_node] = true
    end
    return solution
end

function steepest_local_search_candidate(starting_solution, distance_matrix, costs, neighbors)
    n = length(starting_solution)
    solution = copy(starting_solution)
    objective = calculate_cycle_length(solution, distance_matrix, costs)
    total_nodes = size(distance_matrix, 1)
    
    pos = zeros(Int, total_nodes)
    for (i, node) in enumerate(solution)
        pos[node] = i
    end
    
    while true
        best_delta = -1e-6
        best_move_type = :none
        best_move_args = ()
        
        # 1. Intra-route (2-opt)
        for i in 1:n
            u = solution[i]
            for k in 1:size(neighbors, 2)
                v = neighbors[u, k]
                j = pos[v]
                
                if j != 0
                    a, b = i, j
                    if a > b; a, b = b, a; end
                    
                    if b == a + 1 || (a == 1 && b == n); continue; end
                    
                    delta = intra_two_edges_delta(solution, a, b, distance_matrix, n)
                    if delta < best_delta
                        best_delta = delta
                        best_move_type = :intra
                        best_move_args = (a, b)
                    end
                end
            end
        end
        
        # 2. Inter-route (Swap/Replace)
        for i in 1:n
            prev = solution[i == 1 ? n : i - 1]
            next = solution[i == n ? 1 : i + 1]
            
            # Check neighbors of prev
            for k in 1:size(neighbors, 2)
                v = neighbors[prev, k]
                if pos[v] == 0
                    delta = inter_two_nodes_delta(solution, i, v, distance_matrix, costs, n)
                    if delta < best_delta
                        best_delta = delta
                        best_move_type = :inter
                        best_move_args = (i, v)
                    end
                end
            end
            
            # Check neighbors of next
            for k in 1:size(neighbors, 2)
                v = neighbors[next, k]
                if pos[v] == 0
                    delta = inter_two_nodes_delta(solution, i, v, distance_matrix, costs, n)
                    if delta < best_delta
                        best_delta = delta
                        best_move_type = :inter
                        best_move_args = (i, v)
                    end
                end
            end
        end
        
        if best_move_type == :none
            break
        end
        
        if best_move_type == :intra
            a, b = best_move_args
            apply_intra_two_edges!(solution, a, b)
            for k in (a + 1):b
                pos[solution[k]] = k
            end
        else
            idx, new_node = best_move_args
            old_node = apply_inter_two_nodes!(solution, idx, new_node)
            pos[old_node] = 0
            pos[new_node] = idx
        end
        
        objective += best_delta
    end
    
    return objective, solution
end

function steepest_local_search_candidate_memory(starting_solution, distance_matrix, costs, neighbors)
    n = length(starting_solution)
    solution = copy(starting_solution)
    objective = calculate_cycle_length(solution, distance_matrix, costs)
    total_nodes = size(distance_matrix, 1)
    k_neighbors = size(neighbors, 2)
    
    pos = zeros(Int, total_nodes)
    for (i, node) in enumerate(solution)
        pos[node] = i
    end
    
    deltas = fill(Inf, total_nodes, k_neighbors)
    move_types = fill(:none, total_nodes, k_neighbors)
    move_args = Matrix{Tuple{Int, Int}}(undef, total_nodes, k_neighbors)
    
    # Helper to calculate delta
    function calculate_move(u, k)
        v = neighbors[u, k]
        i = pos[u]
        j = pos[v]
        
        if j != 0 # Intra
            a, b = i, j
            if a > b; a, b = b, a; end
            if b == a + 1 || (a == 1 && b == n)
                return Inf, :none, (0, 0)
            end
            d = intra_two_edges_delta(solution, a, b, distance_matrix, n)
            return d, :intra, (a, b)
        else # Inter
            # We need to check if u is connected to v in the tour?
            # In the original code:
            # Check neighbors of prev: if v is neighbor of prev, try replacing i (current) with v?
            # No, original code:
            # for i in 1:n
            #   prev = ...
            #   for k... v = neighbors[prev, k] ... inter_two_nodes_delta(solution, i, v...)
            
            # This means we iterate over EDGES (prev, i) and check if we can replace i with v where v is neighbor of prev.
            # My memory structure is indexed by u.
            # If I index by u, I can store moves originating from u.
            
            # Original code logic for Inter:
            # Iterate i (current node).
            # Check neighbors of prev (node before i).
            # Check neighbors of next (node after i).
            
            # If I use u as the "anchor", I can store moves.
            # But Inter moves are defined by the node being replaced (i).
            # And the candidate v comes from neighbors of prev or next.
            
            # Let's stick to the original logic for Inter moves.
            # But how to map to memory[u, k]?
            # If u is prev, and v is neighbor of prev. Move replaces i (next of prev).
            # If u is next, and v is neighbor of next. Move replaces i (prev of next).
            
            # So for u, and neighbor v (not in solution):
            # We can consider replacing u_next with v.
            # We can consider replacing u_prev with v.
            
            # So memory[u, k] can store the best Inter move involving u and v?
            # Or we can store multiple moves?
            # Or we can just store the best of the two?
            
            # Let's calculate both and take min.
            
            best_d = Inf
            best_args = (0, 0)
            
            # 1. Replace u_next with v
            u_next_idx = i == n ? 1 : i + 1
            # u is prev of u_next_idx.
            d1 = inter_two_nodes_delta(solution, u_next_idx, v, distance_matrix, costs, n)
            if d1 < best_d
                best_d = d1
                best_args = (u_next_idx, v)
            end
            
            # 2. Replace u_prev with v
            u_prev_idx = i == 1 ? n : i - 1
            # u is next of u_prev_idx.
            d2 = inter_two_nodes_delta(solution, u_prev_idx, v, distance_matrix, costs, n)
            if d2 < best_d
                best_d = d2
                best_args = (u_prev_idx, v)
            end
            
            return best_d, :inter, best_args
        end
    end
    
    # Initialize
    for u in solution
        for k in 1:k_neighbors
            d, t, a = calculate_move(u, k)
            deltas[u, k] = d
            move_types[u, k] = t
            move_args[u, k] = a
        end
    end
    
    affected_mask = falses(total_nodes)
    
    while true
        best_delta = -1e-6
        best_u = -1
        best_k = -1
        
        for u in solution
            for k in 1:k_neighbors
                d = deltas[u, k]
                if d < best_delta
                    best_delta = d
                    best_u = u
                    best_k = k
                end
            end
        end
        
        if best_u == -1
            break
        end
        
        type = move_types[best_u, best_k]
        args = move_args[best_u, best_k]
        
        fill!(affected_mask, false)
        
        if type == :intra
            a, b = args
            apply_intra_two_edges!(solution, a, b)
            
            for k in (a + 1):b
                pos[solution[k]] = k
            end
            
            # Affected: solution[a]...solution[b] (indices a to b)
            # Note: solution has been updated.
            # solution[a] is the anchor.
            # solution[a+1]...solution[b] are the reversed nodes.
            # All of them have changed neighbors/edges.
            
            # Mark a...b
            curr = a
            while true
                affected_mask[solution[curr]] = true
                if curr == b
                    break
                end
                curr = curr == n ? 1 : curr + 1
            end
            
            # Mark a-1
            prev_a = a == 1 ? n : a - 1
            affected_mask[solution[prev_a]] = true
            
            # Mark b+1, b+2
            b1 = b == n ? 1 : b + 1
            b2 = b1 == n ? 1 : b1 + 1
            affected_mask[solution[b1]] = true
            affected_mask[solution[b2]] = true
            
        else # Inter
            idx, new_node = args
            old_node = apply_inter_two_nodes!(solution, idx, new_node)
            pos[old_node] = 0
            pos[new_node] = idx
            
            # Mark idx
            affected_mask[new_node] = true # new_node is at solution[idx]
            
            # Mark idx-1, idx-2
            p1 = idx == 1 ? n : idx - 1
            p2 = p1 == 1 ? n : p1 - 1
            affected_mask[solution[p1]] = true
            affected_mask[solution[p2]] = true
            
            # Mark idx+1, idx+2
            n1 = idx == n ? 1 : idx + 1
            n2 = n1 == n ? 1 : n1 + 1
            affected_mask[solution[n1]] = true
            affected_mask[solution[n2]] = true
            
            # Mark old_node
            affected_mask[old_node] = true
        end
        
        objective += best_delta
        
        # Update memory
        # 1. Rows for affected nodes
        # We iterate all nodes, if affected, update row.
        # If node is not in solution, set to Inf (except old_node which is handled)
        
        # We can iterate over affected nodes directly if we collected them, but mask is fast.
        # But we need to iterate over solution to find affected ones?
        # Or just iterate 1:total_nodes?
        # Iterating solution is faster (n vs N).
        
        # Also need to handle old_node which is NOT in solution anymore.
        # We should set its row to Inf.
        if type == :inter
             # old_node is not in solution, so we don't iterate it in `for u in solution`.
             # But we should clear its entries just in case?
             # Actually, we only read deltas for u in solution.
             # So we don't need to clear old_node's row.
             nothing
        end
        
        # Update rows for nodes in solution that are affected
        for u in solution
            if affected_mask[u]
                for k in 1:k_neighbors
                    d, t, a = calculate_move(u, k)
                    deltas[u, k] = d
                    move_types[u, k] = t
                    move_args[u, k] = a
                end
            end
        end
        
        # 2. Entries where neighbor is affected
        for u in solution
            # If u is affected, we already updated its row.
            if affected_mask[u]
                continue
            end
            
            for k in 1:k_neighbors
                v = neighbors[u, k]
                if affected_mask[v]
                    d, t, a = calculate_move(u, k)
                    deltas[u, k] = d
                    move_types[u, k] = t
                    move_args[u, k] = a
                end
            end
        end
    end
    
    objective = calculate_cycle_length(solution, distance_matrix, costs)
    return objective, solution
end

function double_bridge_perturb(solution; rng=Random.default_rng())
    n = length(solution)
    if n < 8
        return copy(solution)
    end
    
    # Pick 4 split points
    # We need 4 indices such that 1 <= p1 < p2 < p3 < p4 < n
    # Actually we want segments.
    # indices = sort(sample(rng, 1:n-1, 4, replace=false))
    # But sample is in StatsBase. We are using Random.
    # We can use randperm.
    
    indices = sort(randperm(rng, n-1)[1:4])
    p1, p2, p3, p4 = indices[1], indices[2], indices[3], indices[4]
    
    # Segments:
    # 1..p1
    # p1+1..p2
    # p2+1..p3
    # p3+1..p4
    # p4+1..n
    
    # Reorder to S1 - S4 - S3 - S2 - S5
    
    new_solution = Vector{Int}(undef, n)
    
    # S1
    dest = 1
    len = p1
    copyto!(new_solution, dest, solution, 1, len)
    dest += len
    
    # S4 (p3+1..p4)
    len = p4 - p3
    copyto!(new_solution, dest, solution, p3+1, len)
    dest += len
    
    # S3 (p2+1..p3)
    len = p3 - p2
    copyto!(new_solution, dest, solution, p2+1, len)
    dest += len
    
    # S2 (p1+1..p2)
    len = p2 - p1
    copyto!(new_solution, dest, solution, p1+1, len)
    dest += len
    
    # S5 (p4+1..n)
    len = n - p4
    copyto!(new_solution, dest, solution, p4+1, len)
    
    return new_solution
end

function iterated_local_search(distance_matrix, costs, run_time, neighbors; perturb_size=4, rng=Random.default_rng(), initial_solution=nothing, use_memory=false)
    x = isnothing(initial_solution) ? construct_greedy_solution(distance_matrix, costs; rng=rng) : copy(initial_solution)
    
    ls_func = use_memory ? steepest_local_search_candidate_memory : steepest_local_search_candidate

    # Initial Local Search
    objective, x = ls_func(x, distance_matrix, costs, neighbors)
    
    best_objective = objective
    best_solution = copy(x)
    
    start_time = time()
    i = 0

    while (time() - start_time) <= run_time
        i += 1
        # Use Double Bridge Perturbation
        y = double_bridge_perturb(x; rng=rng)
        # Add small extra noise to break structures preserved by double bridge
        if rand(rng) < 0.5
            y = perturb(y, distance_matrix, 1; rng=rng)
        end
        
        objective_y, y = ls_func(y, distance_matrix, costs, neighbors)
        
        if objective_y < objective
            objective, x = objective_y, y
            if objective < best_objective
                best_objective = objective
                best_solution = copy(x)
            end
        elseif rand(rng) < 0.01 # Small chance to accept worse solution
             objective, x = objective_y, y
        end
    end

    return best_objective, best_solution, i
end

function multiple_start_local_search(distance_matrix, costs; rng=Random.default_rng())
    objective = Inf
    solution = Vector{Int}()
    for _ in 1:200
        start_solution = random_start(distance_matrix; rng=rng)
        new_objective, new_solution = steepest_local_search(start_solution, distance_matrix, costs)

        if new_objective < objective
            objective = new_objective
            solution = new_solution
        end
    end

    return objective, solution
end

function steepest_local_search(starting_solution, distance_matrix, costs)
    n = length(starting_solution)
    solution = copy(starting_solution)
    objective = calculate_cycle_length(solution, distance_matrix, costs)
    total_nodes = size(distance_matrix, 1)

    intra_moves = precompute_intra_moves(n)
    available_nodes, node_positions = initialize_available_nodes(solution, total_nodes)

    while true
        best_delta = 0.0
        best_kind = 0  # 0: none, 1: intra, 2: inter
        best_a = 0
        best_b = 0

        @inbounds for (a, b) in intra_moves
            delta = intra_two_edges_delta(solution, a, b, distance_matrix, n)
            if delta < best_delta
                best_delta = delta
                best_kind = 1
                best_a = a
                best_b = b
            end
        end

        @inbounds for a in 1:n
            for idx in eachindex(available_nodes)
                b = available_nodes[idx]
                delta = inter_two_nodes_delta(solution, a, b, distance_matrix, costs, n)
                if delta < best_delta
                    best_delta = delta
                    best_kind = 2
                    best_a = a
                    best_b = b
                end
            end
        end

        if best_kind == 0
            break
        end

        if best_kind == 1
            apply_intra_two_edges!(solution, best_a, best_b)
        else
            old_node = apply_inter_two_nodes!(solution, best_a, best_b)
            remove_available!(available_nodes, node_positions, best_b)
            add_available!(available_nodes, node_positions, old_node)
        end

        objective += best_delta
    end

    return objective, solution
end

function precompute_intra_moves(n)
    moves = Vector{NTuple{2,Int}}()
    sizehint!(moves, max(0, n * (n - 3) ÷ 2))
    for i in 1:n
        for j in i+2:n
            if abs(i - j) != n - 1
                push!(moves, (i, j))
            end
        end
    end
    return moves
end

function initialize_available_nodes(solution, total_nodes)
    available_nodes = Vector{Int}()
    sizehint!(available_nodes, total_nodes - length(solution))
    node_positions = zeros(Int, total_nodes)
    in_solution = falses(total_nodes)
    @inbounds for node in solution
        in_solution[node] = true
    end

    for node in 1:total_nodes
        if !in_solution[node]
            push!(available_nodes, node)
            node_positions[node] = length(available_nodes)
        end
    end

    return available_nodes, node_positions
end

@inline function remove_available!(available_nodes, node_positions, node)
    idx = node_positions[node]
    if idx == 0
        return
    end
    last_idx = length(available_nodes)
    if idx != last_idx
        swap_node = available_nodes[last_idx]
        available_nodes[idx] = swap_node
        node_positions[swap_node] = idx
    end
    pop!(available_nodes)
    node_positions[node] = 0
end

@inline function add_available!(available_nodes, node_positions, node)
    push!(available_nodes, node)
    node_positions[node] = length(available_nodes)
end