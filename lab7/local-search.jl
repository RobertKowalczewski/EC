using Random, ProgressMeter, Statistics

include("helpers.jl")
include("moves.jl")


fmt(x) = round(x; digits=5)

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