using CSV
using DataFrames
using Random

include("moves.jl")

struct TSPInstance
    distance::Matrix{Int}
    costs::Vector{Int}
    required::Int
end

Base.length(instance::TSPInstance) = size(instance.distance, 1)

mutable struct TSPSolution
    tour::Vector{Int}
    selected::BitVector
    tour_length::Int
    node_cost_sum::Int
    objective::Int
end

mutable struct SearchStats
    iterations::Int
    evaluated_moves::Int
    applied_moves::Int
end

SearchStats() = SearchStats(0, 0, 0)

function load_instance(path::AbstractString)
    data = CSV.read(path, DataFrame; header=["x", "y", "w"])
    n = nrow(data)
    distance = Matrix{Int}(undef, n, n)
    xs = Float64.(data[:, :x])
    ys = Float64.(data[:, :y])
    for i in 1:n
        distance[i, i] = 0
        for j in i + 1:n
            dx = xs[i] - xs[j]
            dy = ys[i] - ys[j]
            dist = round(Int, sqrt(dx * dx + dy * dy))
            distance[i, j] = dist
            distance[j, i] = dist
        end
    end
    costs = Int.(data[:, :w])
    required = ceil(Int, n / 2)
    return TSPInstance(distance, costs, required)
end

function cycle_distance(tour::Vector{Int}, distance::Matrix{Int})
    n = length(tour)
    n == 0 && return 0
    total = 0
    for idx in 1:n
        from = tour[idx]
        to = tour[mod1(idx + 1, n)]
        total += distance[from, to]
    end
    return total
end

objective(solution::TSPSolution) = solution.objective

function build_solution(tour::Vector{Int}, instance::TSPInstance)
    n = length(instance)
    selected = falses(n)
    for node in tour
        selected[node] = true
    end
    tour_length = cycle_distance(tour, instance.distance)
    node_cost_sum = sum(instance.costs[tour])
    return TSPSolution(copy(tour), selected, tour_length, node_cost_sum, tour_length + node_cost_sum)
end

function random_initial_solution(instance::TSPInstance, rng::AbstractRNG)
    nodes = collect(1:length(instance))
    shuffle!(rng, nodes)
    tour = nodes[1:instance.required]
    shuffle!(rng, tour)
    return build_solution(tour, instance)
end

function two_regret_greedy_initial_solution(instance::TSPInstance, start::Int)
    n = length(instance)
    required = instance.required
    distance = instance.distance
    costs = instance.costs

    tour = Int[start]
    visited = falses(n)
    visited[start] = true

    if required == 1
        return build_solution(tour, instance)
    end

    best_second = 0
    best_score = typemax(Int)
    for node in 1:n
        if visited[node]
            continue
        end
        score = distance[start, node] + distance[node, start] + costs[node]
        if score < best_score
            best_score = score
            best_second = node
        end
    end
    if best_second == 0
        error("Unable to select second node for greedy initialization")
    end

    push!(tour, best_second)
    visited[best_second] = true

    while length(tour) < required
        best_candidate = 0
        best_position = 0
        best_cost = typemax(Int)
        best_score = typemin(Int)

        for node in 1:n
            if visited[node]
                continue
            end
            best_insert_cost = typemax(Int)
            second_insert_cost = typemax(Int)
            best_insert_pos = 0
            for idx in 1:length(tour)
                from = tour[idx]
                to = tour[mod1(idx + 1, length(tour))]
                insert_cost = distance[from, node] + distance[node, to] - distance[from, to] + costs[node]
                if insert_cost < best_insert_cost
                    second_insert_cost = best_insert_cost
                    best_insert_cost = insert_cost
                    best_insert_pos = idx + 1
                elseif insert_cost < second_insert_cost
                    second_insert_cost = insert_cost
                end
            end
            if second_insert_cost == typemax(Int)
                second_insert_cost = best_insert_cost
            end
            regret = second_insert_cost - best_insert_cost
            score = regret
            if score > best_score || (score == best_score && best_insert_cost < best_cost)
                best_score = score
                best_cost = best_insert_cost
                best_candidate = node
                best_position = best_insert_pos
            end
        end

        if best_candidate == 0 || best_position == 0
            error("Failed to determine insertion during greedy initialization")
        end

        insert!(tour, best_position, best_candidate)
        visited[best_candidate] = true
    end

    return build_solution(tour, instance)
end

function copy_solution(solution::TSPSolution)
    return TSPSolution(copy(solution.tour), copy(solution.selected), solution.tour_length, solution.node_cost_sum, solution.objective)
end

function node_after_swap(tour::Vector{Int}, i::Int, j::Int, idx::Int)
    if idx == i
        return tour[j]
    elseif idx == j
        return tour[i]
    else
        return tour[idx]
    end
end

function delta_intra_node(solution::TSPSolution, instance::TSPInstance, move::IntraNodeMove)
    n = length(solution.tour)
    n <= 1 && return (0, 0)
    i, j = move.i, move.j
    if i == j
        return (0, 0)
    end
    if j < i
        i, j = j, i
    end

    affected = Int[mod1(i - 1, n), i, mod1(j - 1, n), j]
    sort!(affected)
    unique!(affected)
    distance = instance.distance
    tour = solution.tour

    delta = 0
    for idx in affected
        from_old = tour[idx]
        to_index = mod1(idx + 1, n)
        to_old = tour[to_index]
        new_from = node_after_swap(tour, move.i, move.j, idx)
        new_to = node_after_swap(tour, move.i, move.j, to_index)
        delta += distance[new_from, new_to] - distance[from_old, to_old]
    end

    return (delta, 0)
end

function delta_intra_edge(solution::TSPSolution, instance::TSPInstance, move::IntraEdgeMove)
    n = length(solution.tour)
    n < 2 && return (0, 0)
    i, j = move.i, move.j
    if i == j
        return (0, 0)
    end
    if j < i
        i, j = j, i
    end
    if i == 1 && j == n
        return (0, 0)
    end

    tour = solution.tour
    distance = instance.distance
    prev_i = tour[mod1(i - 1, n)]
    node_i = tour[i]
    node_j = tour[j]
    next_j = tour[mod1(j + 1, n)]

    before = distance[prev_i, node_i] + distance[node_j, next_j]
    after = distance[prev_i, node_j] + distance[node_i, next_j]
    return (after - before, 0)
end

function delta_inter_node(solution::TSPSolution, instance::TSPInstance, move::InterNodeMove)
    n = length(solution.tour)
    n == 0 && return (0, 0)
    pos = move.position
    new_node = move.node
    old_node = solution.tour[pos]
    if new_node == old_node
        return (0, 0)
    end
    distance = instance.distance
    prev = solution.tour[mod1(pos - 1, n)]
    nxt = solution.tour[mod1(pos + 1, n)]

    before = distance[prev, old_node] + distance[old_node, nxt]
    after = distance[prev, new_node] + distance[new_node, nxt]
    delta_distance = after - before
    delta_cost = instance.costs[new_node] - instance.costs[old_node]
    return (delta_distance, delta_cost)
end

function calculate_delta(solution::TSPSolution, instance::TSPInstance, move::AbstractMove)
    if move isa IntraNodeMove
        return delta_intra_node(solution, instance, move)
    elseif move isa IntraEdgeMove
        return delta_intra_edge(solution, instance, move)
    elseif move isa InterNodeMove
        return delta_inter_node(solution, instance, move)
    else
        error("Unsupported move type")
    end
end

function apply_move!(solution::TSPSolution, instance::TSPInstance, move::IntraNodeMove, delta::Tuple{Int, Int})
    solution.tour[move.i], solution.tour[move.j] = solution.tour[move.j], solution.tour[move.i]
    solution.tour_length += delta[1]
    solution.node_cost_sum += delta[2]
    solution.objective += delta[1] + delta[2]
end

function apply_move!(solution::TSPSolution, instance::TSPInstance, move::IntraEdgeMove, delta::Tuple{Int, Int})
    i, j = move.i, move.j
    if j < i
        i, j = j, i
    end
    reverse!(solution.tour, i, j)
    solution.tour_length += delta[1]
    solution.node_cost_sum += delta[2]
    solution.objective += delta[1] + delta[2]
end



function apply_move!(solution::TSPSolution, instance::TSPInstance, move::AbstractMove, delta::Tuple{Int, Int})
    if n < 2
        return AbstractMove[]
    end
    if move isa IntraNodeMove
        apply_move!(solution, instance, move::IntraNodeMove, delta)
    elseif move isa IntraEdgeMove
        apply_move!(solution, instance, move::IntraEdgeMove, delta)
    elseif move isa InterNodeMove
        apply_move!(solution, instance, move::InterNodeMove, delta)
    end
end
        if n < 3
            return moves
        end

function generate_intra_moves(solution::TSPSolution, intra_type::Symbol)
    n = length(solution.tour)
    moves = AbstractMove[]
    if intra_type == :node
        for i in 1:n-1
            for j in i + 1:n
                push!(moves, IntraNodeMove(i, j))
            end
        end
    elseif intra_type == :edge
        for i in 1:n-1
            for j in i + 1:n
                if i == 1 && j == n
                    continue
                end
                push!(moves, IntraEdgeMove(i, j))
            end
        end
    else
        error("Unsupported intra-route move type")
    end
    return moves
end

function generate_inter_moves(solution::TSPSolution)
    nsel = length(solution.tour)
    moves = InterNodeMove[]
    if nsel == 0
        return moves
    end
    unselected = findall(!, solution.selected)
    for pos in 1:nsel
        for node in unselected
            push!(moves, InterNodeMove(pos, node))
        end
    end
    return moves
end

function all_moves(solution::TSPSolution, intra_type::Symbol)
    moves = AbstractMove[]
    append!(moves, generate_intra_moves(solution, intra_type))
    append!(moves, generate_inter_moves(solution))
    return moves
end

function recompute_objective(solution::TSPSolution, instance::TSPInstance)
    tour_length = cycle_distance(solution.tour, instance.distance)
    node_cost_sum = sum(instance.costs[solution.tour])
    return tour_length + node_cost_sum
end