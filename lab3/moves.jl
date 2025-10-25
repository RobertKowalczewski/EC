function intra_two_nodes_exchange(solution, objective, a, b, distance_matrix)
    n = length(solution)

    if a > b
        a, b = b, a
    end

    prev_a = a == 1 ? n : a - 1
    next_a = a == n ? 1 : a + 1
    prev_b = b == 1 ? n : b - 1
    next_b = b == n ? 1 : b + 1

    if (a == 1 && b == n)
        old_edges = distance_matrix[solution[b], solution[a]] + distance_matrix[solution[prev_b], solution[b]] + distance_matrix[solution[a], solution[next_a]]
        new_edges = distance_matrix[solution[prev_b], solution[a]] + distance_matrix[solution[a], solution[b]] + distance_matrix[solution[b], solution[next_a]]

    elseif abs(a-b) == 1
        old_edges = distance_matrix[solution[prev_a], solution[a]] + distance_matrix[solution[a], solution[b]] + distance_matrix[solution[b], solution[next_b]]
        new_edges = distance_matrix[solution[prev_a], solution[b]] + distance_matrix[solution[b], solution[a]] + distance_matrix[solution[a], solution[next_b]]
    else
        old_edges = distance_matrix[solution[prev_a], solution[a]] + distance_matrix[solution[a], solution[next_a]] +
                    distance_matrix[solution[prev_b], solution[b]] + distance_matrix[solution[b], solution[next_b]]
        new_edges = distance_matrix[solution[prev_a], solution[b]] + distance_matrix[solution[b], solution[next_a]] +
                     distance_matrix[solution[prev_b], solution[a]] +  distance_matrix[solution[a], solution[next_b]]
    end

    delta = new_edges - old_edges  # cost change in path only
    new_objective = objective + delta  # node costs unchanged

    new_solution = copy(solution)
    new_solution[a], new_solution[b] = new_solution[b], new_solution[a]

    return new_objective, new_solution

end

function intra_two_edges_exchange(solution, objective, a, b, distance_matrix)

end

function inter_two_nodes_exchange(solution, objective, a, b, distance_matrix, costs)

end