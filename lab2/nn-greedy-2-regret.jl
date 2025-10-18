function nn_greedy_2_regret(distance_matrix, costs, weights, start=1)
    n = size(distance_matrix)[1]
    cycle_length = ceil(Int, n / 2)
    solution = [start]
    visited = Set([start])

    # Add the second node: choose the nearest neighbor to start
    best_distance, best_node = Inf, nothing
    for i in 1:n
        if i ∉ visited
            d = distance_matrix[start, i] + costs[i]
            if d < best_distance
                best_distance = d
                best_node = i
            end
        end
    end
    push!(solution, best_node)
    push!(visited, best_node)

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

                # insert before the first node
                insertion_cost = distance_matrix[i, solution[1]] + costs[i]
                best_cost = insertion_cost
                best_pos = 1

                # insert after the last node
                insertion_cost = distance_matrix[solution[end], i] + costs[i]
                if insertion_cost < best_cost
                    second_best_cost = best_cost

                    best_cost = insertion_cost
                    best_pos = length(solution) + 1
                else
                    second_best_cost = insertion_cost
                end

                # Find insertion cost for each possible position in the current solution (excluding before first and after last node)
                for pos in 1:(length(solution)-1)
                    prev = solution[pos]
                    next = solution[pos+1]
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