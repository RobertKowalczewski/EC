using Random, ProgressMeter, Statistics

include("moves.jl")
include("helpers.jl")

const Edge = Tuple{Int,Int}
const EdgeInfo = NamedTuple{(:orientation, :start_idx, :end_idx),Tuple{Symbol,Int,Int}}

struct SavedMove
    delta::Float64
    move_kind::Symbol
    old_edges::Vector{Edge}
    new_edges::Vector{Edge}
    orientations::Vector{Symbol}
end

function heap_push!(heap::Vector{SavedMove}, move::SavedMove)
    push!(heap, move)
    idx = length(heap)
    while idx > 1
        parent = idx >>> 1
        if heap[parent].delta <= heap[idx].delta
            break
        end
        heap[parent], heap[idx] = heap[idx], heap[parent]
        idx = parent
    end
end

function heap_pop!(heap::Vector{SavedMove})
    isempty(heap) && return nothing
    top = heap[1]
    last = pop!(heap)
    if isempty(heap)
        return top
    end
    heap[1] = last
    idx = 1
    len = length(heap)
    while true
        left = idx << 1
        right = left + 1
        smallest = idx
        if left <= len && heap[left].delta < heap[smallest].delta
            smallest = left
        end
        if right <= len && heap[right].delta < heap[smallest].delta
            smallest = right
        end
        smallest == idx && break
        heap[idx], heap[smallest] = heap[smallest], heap[idx]
        idx = smallest
    end
    return top
end


function test_algorithms(distance_matrix, costs, data_path, data_name)
    for start_type in [random_start]
        for local_search_func in [steepest_local_search]
            for intra_function in [intra_two_edges_exchange]

                println(start_type, " ", local_search_func, " ", intra_function)

                objectives = []
                solutions = []
                starts = []
                start_objectives = []
                times = []

                fmt(x) = round(x; digits=5)

                @showprogress 0.1 "runs" for run in 1:200
                    starting_solution = start_type(distance_matrix, costs, [0.5, 0.5], run)
                    t = @elapsed begin
                        objective, solution = local_search_func(starting_solution, distance_matrix, costs, intra_function, inter_two_nodes_exchange)
                    end
                    push!(times, t)
                    push!(objectives, objective)
                    push!(solutions, solution)

                    push!(starts, starting_solution)
                    push!(start_objectives, calculate_cycle_length(starting_solution, distance_matrix, costs))

                    if run == 1
                        println("solution from first node:")
                        println([x - 1 for x in solution])
                    end
                end

                min_s = argmin(objectives)
                max_s = argmax(objectives)
                avg_s = mean(objectives)

                min_s_start = minimum(start_objectives)
                max_s_start = maximum(start_objectives)
                avg_s_start = mean(start_objectives)

                avg_time = mean(times)
                min_time = minimum(times)
                max_time = maximum(times)


                println("After Local Search | Starting Solution")
                println("Obj: $(fmt(avg_s)) ($(fmt(objectives[min_s])) - $(fmt(objectives[max_s]))) | Obj: $(fmt(avg_s_start)) ($(fmt(min_s_start)) - $(fmt(max_s_start)))")
                println("Time[s]: $(avg_time) ($(min_time) - $(max_time))")

                plot_best_solution(objectives, solutions,
                    data_path,
                    "$(start_type)_$(local_search_func)_$(intra_function)",
                    "lab5v2/$(data_name)/$(start_type)_$(local_search_func)_$(intra_function).png")
            end
        end
    end
end

edge_step_forward(idx, n) = idx == n ? 1 : idx + 1

function build_position_map(solution)
    mapping = Dict{Int,Int}()
    for (idx, node) in enumerate(solution)
        mapping[node] = idx
    end
    return mapping
end

function locate_edge(solution, edge::Edge, position_map)
    n = length(solution)
    a, b = edge

    idx_a = get(position_map, a, nothing)
    if !isnothing(idx_a)
        next_idx = edge_step_forward(idx_a, n)
        solution[next_idx] == b && return (:forward, idx_a, next_idx)
    end

    idx_b = get(position_map, b, nothing)
    if !isnothing(idx_b)
        next_idx = edge_step_forward(idx_b, n)
        solution[next_idx] == a && return (:reverse, idx_b, next_idx)
    end

    return (:missing, nothing, nothing)
end

function capture_edge_orientations(solution, old_edges, position_map)
    orientations = Symbol[]
    for edge in old_edges
        orientation, _, _ = locate_edge(solution, edge, position_map)
        orientation == :missing && return nothing
        push!(orientations, orientation)
    end
    return orientations
end

function collect_circular_indices(start_idx, end_idx, n)
    indices = Int[]
    idx = start_idx
    push!(indices, idx)
    while idx != end_idx
        idx = edge_step_forward(idx, n)
        if length(indices) > n
            return nothing
        end
        push!(indices, idx)
    end
    return indices
end

function apply_intra_move(solution, edge_infos)
    length(edge_infos) == 2 || return nothing
    first_edge, second_edge = edge_infos
    n = length(solution)
    segment_indices = collect_circular_indices(first_edge.end_idx, second_edge.start_idx, n)
    isnothing(segment_indices) && return nothing

    new_solution = copy(solution)
    segment_nodes = solution[segment_indices]
    reverse!(segment_nodes)

    for (idx, node) in zip(segment_indices, segment_nodes)
        new_solution[idx] = node
    end

    return new_solution
end

function find_repeated_node(edge_list)
    counts = Dict{Int,Int}()
    for (u, v) in edge_list
        counts[u] = get(counts, u, 0) + 1
        counts[v] = get(counts, v, 0) + 1
    end
    for (node, count) in counts
        count > 1 && return node
    end
    return nothing
end

function apply_inter_move(solution, move::SavedMove, position_map)
    removed_node = find_repeated_node(move.old_edges)
    new_node = find_repeated_node(move.new_edges)

    if removed_node === nothing || new_node === nothing
        return nothing
    end

    idx = get(position_map, removed_node, nothing)
    isnothing(idx) && return nothing
    haskey(position_map, new_node) && return nothing

    new_solution = copy(solution)
    new_solution[idx] = new_node
    return new_solution
end

function check_if_applicable(solution, position_map, move::SavedMove)
    length(move.old_edges) == length(move.orientations) || return (:remove, nothing)
    edge_infos = EdgeInfo[]
    for edge in move.old_edges
        orientation, start_idx, end_idx = locate_edge(solution, edge, position_map)
        if orientation == :missing
            return (:remove, nothing)
        end
        push!(edge_infos, (orientation=orientation, start_idx=start_idx, end_idx=end_idx))
    end

    matches = true
    reversed = true
    for (idx, info) in enumerate(edge_infos)
        stored_orientation = move.orientations[idx]
        matches &= info.orientation == stored_orientation
        reversed &= info.orientation != stored_orientation
    end

    if !(matches || reversed)
        return (:skip, nothing)
    end

    new_solution = move.move_kind == :intra ?
                   apply_intra_move(solution, edge_infos) :
                   apply_inter_move(solution, move, position_map)

    new_solution === nothing && return (:remove, nothing)
    return (:apply, new_solution)
end

function pick_cached_move!(LM, solution, position_map)
    skipped = SavedMove[]

    while !isempty(LM)
        move = heap_pop!(LM)
        status, candidate = check_if_applicable(solution, position_map, move)
        if status == :remove
            continue
        elseif status == :skip
            push!(skipped, move)
        else
            for frozen in skipped
                heap_push!(LM, frozen)
            end
            return move, candidate
        end
    end

    for frozen in skipped
        heap_push!(LM, frozen)
    end

    return nothing, nothing
end

function evaluate_new_moves!(LM, solution, objective, intra_move, distance_matrix, costs,
    intra_function, inter_function, position_map)
    inserted = 0

    for (a, b) in intra_move
        new_objective, old_edges, new_edges = intra_function(solution, objective, a, b, distance_matrix)
        isempty(old_edges) && continue
        delta = new_objective - objective
        delta >= 0 && continue

        orientations = capture_edge_orientations(solution, old_edges, position_map)
        orientations === nothing && continue

        heap_push!(LM, SavedMove(delta, :intra, old_edges, new_edges, orientations))
        inserted += 1
    end

    available_nodes = setdiff(1:size(distance_matrix, 1), solution)
    sol_len = length(solution)
    for a in 1:sol_len
        for b in available_nodes
            new_objective, old_edges, new_edges = inter_function(solution, objective, a, b, distance_matrix, costs)
            isempty(old_edges) && continue
            delta = new_objective - objective
            delta >= 0 && continue

            orientations = capture_edge_orientations(solution, old_edges, position_map)
            orientations === nothing && continue

            heap_push!(LM, SavedMove(delta, :inter, old_edges, new_edges, orientations))
            inserted += 1
        end
    end

    return inserted
end

function steepest_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)

    intra_move = [(i, j) for i in 1:n for j in i+2:n if abs(i - j) != n - 1]
    LM = SavedMove[]
    position_map = build_position_map(solution)

    while true
        selected_move, candidate_solution = pick_cached_move!(LM, solution, position_map)
        if selected_move !== nothing
            solution = candidate_solution
            objective += selected_move.delta
            position_map = build_position_map(solution)
            continue
        end

        inserted = evaluate_new_moves!(LM, solution, objective, intra_move, distance_matrix, costs,
            intra_function, inter_function, position_map)

        inserted == 0 && break
    end

    return objective, solution
end