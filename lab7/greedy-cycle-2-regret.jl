function greedy_cycle_2_regret(distance_matrix, costs, weights, incomplete_solution)
    n = size(distance_matrix)[1]
    cycle_length = ceil(Int, n / 2)
    solution = copy(incomplete_solution)
    visited = Set(solution)

    # # Add the second node: choose the nearest neighbor to start
    # best_distance, best_node = Inf, nothing
    # for i in 1:n
    #     if i ∉ visited
    #         d = distance_matrix[start, i] + costs[i]
    #         if d < best_distance
    #             best_distance = d
    #             best_node = i
    #         end
    #     end
    # end
    # push!(solution, best_node)
    # push!(visited, best_node)

    while length(solution) < cycle_length
        best_score = -Inf
        best_node = nothing
        best_insertion_pos = nothing

        for i in 1:n
            if i ∉ visited
                # Track positions and costs
                best_cost = Inf
                second_best_cost = Inf
                best_pos = nothing

                # Find insertion cost for each possible position in the current solution
                for pos in eachindex(solution)
                    prev = solution[pos]
                    next = solution[mod1(pos + 1, length(solution))]
                    insertion_cost = distance_matrix[prev, i] + costs[i] + distance_matrix[i, next] - distance_matrix[prev, next]

                    if insertion_cost < best_cost
                        second_best_cost = best_cost
                        best_cost = insertion_cost
                        best_pos = pos + 1
                    elseif insertion_cost < second_best_cost
                        second_best_cost = insertion_cost
                    end
                end

                regret = second_best_cost - best_cost
                score = weights[1] * regret - weights[2] * best_cost # if weights=[1,0], weighted 2-regret heuristic becomes identical to the pure greedy 2-regret heuristic

                if score > best_score
                    best_score = score
                    best_node = i
                    best_insertion_pos = best_pos
                end
            end
        end

        insert!(solution, best_insertion_pos, best_node)
        push!(visited, best_node)
    end

    return solution
end