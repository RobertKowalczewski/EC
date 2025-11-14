
function intra_two_edges_exchange(solution, objective, a, b, distance_matrix)
    n = length(solution)
    if a == b || abs(a - b) == 1 || abs(a - b) == n - 1
        println("INVALID NODES")
        return objective, solution
    end
    a1 = a
    a2 = a == n ? 1 : a + 1
    b1 = b
    b2 = b == n ? 1 : b + 1

    old_edges = distance_matrix[solution[a1], solution[a2]] + distance_matrix[solution[b1], solution[b2]]
    new_edges = distance_matrix[solution[a1], solution[b1]] + distance_matrix[solution[a2], solution[b2]]

    delta = new_edges - old_edges
    new_objective = objective + delta  # node costs unchanged

    old_edges = [(solution[a1], solution[a2]), (solution[b1], solution[b2])]
    new_edges = [(solution[a1], solution[b1]), (solution[a2], solution[b2])]

    return new_objective, old_edges, new_edges
end

function inter_two_nodes_exchange(solution, objective, a, b, distance_matrix, costs)
    n = length(solution)

    prev_a = a == 1 ? n : a - 1
    next_a = a == n ? 1 : a + 1

    # Remove node i and insert b in its place
    old_edges = distance_matrix[solution[prev_a], solution[a]] + distance_matrix[solution[a], solution[next_a]]
    new_edges = distance_matrix[solution[prev_a], b] + distance_matrix[b, solution[next_a]]


    # Node cost change
    delta_cost = costs[b] - costs[solution[a]]

    delta = (new_edges - old_edges) + delta_cost
    new_objective = objective + delta

    old_edges = [(solution[prev_a], solution[a]), (solution[a], solution[next_a])]
    new_edges = [(solution[prev_a], b), (b, solution[next_a])]

    return new_objective, old_edges, new_edges
end