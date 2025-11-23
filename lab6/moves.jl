
@inline function intra_two_edges_delta(solution, a, b, distance_matrix, n)
    a2 = a == n ? 1 : a + 1
    b2 = b == n ? 1 : b + 1

    old_edges = distance_matrix[solution[a], solution[a2]] + distance_matrix[solution[b], solution[b2]]
    new_edges = distance_matrix[solution[a], solution[b]] + distance_matrix[solution[a2], solution[b2]]
    return new_edges - old_edges
end

@inline function inter_two_nodes_delta(solution, a, new_node, distance_matrix, costs, n)
    prev_a = a == 1 ? n : a - 1
    next_a = a == n ? 1 : a + 1
    current_node = solution[a]

    old_edges = distance_matrix[solution[prev_a], current_node] + distance_matrix[current_node, solution[next_a]]
    new_edges = distance_matrix[solution[prev_a], new_node] + distance_matrix[new_node, solution[next_a]]
    return (new_edges - old_edges) + (costs[new_node] - costs[current_node])
end

@inline function apply_intra_two_edges!(solution, a, b)
    left = a + 1
    right = b
    while left < right
        solution[left], solution[right] = solution[right], solution[left]
        left += 1
        right -= 1
    end
end

@inline function apply_inter_two_nodes!(solution, a, new_node)
    old_node = solution[a]
    solution[a] = new_node
    return old_node
end