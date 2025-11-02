function precompute_candidate_moves(distance_matrix, costs, n_candidates)
    n = size(distance_matrix, 1)
    candidates = [[] for k in 1:n]
    
    for i in 1:n-1
        for j in i+1:n
            push!(candidates[i], (j, distance_matrix[i, j] + costs[j]))
            push!(candidates[j], (i, distance_matrix[j, i] + costs[i]))
        end
    end
    
    for k in keys(candidates)
        sort!(candidates[k], by=x -> x[2])
    end
    
    return [candidate[1:n_candidates] for candidate in candidates]
end

function steepest_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)
    improvement = true

    if intra_function == intra_two_nodes_exchange
        intra_move = [(i, j) for i in 1:n for j in i+1:n]
    elseif intra_function == intra_two_edges_exchange
        intra_move = [(i, j) for i in 1:n for j in i+2:n if abs(i - j) != n - 1]
    else
        println("invalid intra function")
        return
    end

    while improvement
        improvement = false

        available_nodes = setdiff(1:size(distance_matrix, 1), solution)
        inter_move = [(i, j) for i in 1:n for j in available_nodes]

        moves = [(m, :intra) for m in intra_move]
        append!(moves, [(m, :inter) for m in inter_move])

        best_objective, best_solution = Inf, nothing
        for ((a,b), move_kind) in moves
            if move_kind == :intra
                new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
            else
                new_objective, new_solution = inter_function(solution, objective, a, b, distance_matrix, costs)
            end
            if new_objective < objective && new_objective < best_objective
                best_solution = new_solution
                best_objective = new_objective
                improvement = true
            end
        end

        if improvement
            solution, objective = best_solution, best_objective
        end
    end

    return objective, solution
end


# Directed candidate lists by s(u,v)=d(u,v)+c(v)
function precompute_candidates(distance_matrix, costs; n_candidates::Int=10)
    n = size(distance_matrix, 1)

    cand_out = Vector{Vector{Int}}(undef, n)
    for u in 1:n
        neigh = [(v, distance_matrix[u, v] + costs[v]) for v in 1:n if v != u]
        sort!(neigh, by = x -> x[2])
        k = min(n_candidates, length(neigh))
        cand_out[u] = [neigh[i][1] for i in 1:k]
    end

    cand_in = [Int[] for _ in 1:n]
    for u in 1:n
        for v in cand_out[u]
            push!(cand_in[v], u)
        end
    end

    return cand_out, cand_in
end

# Candidate-restricted steepest local search
function steepest_local_search_candidates(
    starting_solution::Vector{Int},
    distance_matrix,
    costs,
    intra_function,
    inter_function;
    n_candidates::Int=10
)
    n = length(starting_solution)
    solution = copy(starting_solution)
    objective = calculate_cycle_length(solution, distance_matrix, costs)

    cand_out, cand_in = precompute_candidates(distance_matrix, costs; n_candidates=n_candidates)

    # Build node -> position map
    pos_of(sol) = Dict(node => idx for (idx, node) in enumerate(sol))

    improvement = true
    while improvement
        improvement = false
        best_objective = objective
        best_solution = solution

        pos = pos_of(solution)
        solution_set = Set(solution)

        # Inter moves: replace solution[a] with b (unselected) if b introduces a candidate edge
        for a in 1:n
            prev_a = a == 1 ? n : a - 1
            next_a = a == n ? 1 : a + 1
            u = solution[prev_a]
            v = solution[next_a]

            # At least one candidate edge: (u->b) or (b->v)
            candidate_bs = union(cand_out[u], cand_in[v])
            for b in candidate_bs
                if b ∈ solution_set
                    continue
                end
                new_objective, new_solution = inter_function(solution, objective, a, b, distance_matrix, costs)
                if new_objective < best_objective
                    best_objective = new_objective
                    best_solution = new_solution
                    improvement = true
                end
            end
        end

        # Intra moves
        if intra_function == intra_two_edges_exchange
            # 2-opt: (u,v),(w,x) -> (u,w),(v,x); candidate if w ∈ cand_out[u] or x ∈ cand_out[v]
            seen = Set{Tuple{Int,Int}}()
            for i in 1:n
                i2 = (i == n ? 1 : i + 1)
                u = solution[i]
                v = solution[i2]

                # Candidate via (u -> w)
                for w in cand_out[u]
                    if w ∉ solution_set; continue; end
                    j = pos[w]
                    a, b = min(i, j), max(i, j)
                    # Skip adjacent and wrap-adjacent
                    if abs(a - b) == 1 || abs(a - b) == n - 1
                        continue
                    end
                    key = (a, b)
                    if key ∈ seen; continue; end
                    push!(seen, key)

                    new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
                    if new_objective < best_objective
                        best_objective = new_objective
                        best_solution = new_solution
                        improvement = true
                    end
                end

                # Candidate via (v -> x)
                for x in cand_out[v]
                    if x ∉ solution_set; continue; end
                    j2 = pos[x]
                    j1 = (j2 == 1 ? n : j2 - 1)  # edge (j1 -> j2)
                    a, b = min(i, j1), max(i, j1)
                    if a == b || abs(a - b) == 1 || abs(a - b) == n - 1
                        continue
                    end
                    key = (a, b)
                    if key ∈ seen; continue; end
                    push!(seen, key)

                    new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
                    if new_objective < best_objective
                        best_objective = new_objective
                        best_solution = new_solution
                        improvement = true
                    end
                end
            end
        elseif intra_function == intra_two_nodes_exchange
            # Swap positions a,b; candidate if any newly created arc is candidate
            seen = Set{Tuple{Int,Int}}()
            for a in 1:n
                prev_a = a == 1 ? n : a - 1
                next_a = a == n ? 1 : a + 1
                u = solution[prev_a]
                nxt = solution[next_a]

                Bs = union([b for b in cand_out[u] if b ∈ solution_set],
                           [b for b in cand_in[nxt] if b ∈ solution_set])

                for node_b in Bs
                    b = pos[node_b]
                    if b == a; continue; end
                    key = (min(a, b), max(a, b))
                    if key ∈ seen; continue; end
                    push!(seen, key)

                    new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
                    if new_objective < best_objective
                        best_objective = new_objective
                        best_solution = new_solution
                        improvement = true
                    end
                end
            end
        else
            println("invalid intra function")
            return
        end

        if improvement
            solution = best_solution
            objective = best_objective
        end
    end

    return objective, solution
end