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