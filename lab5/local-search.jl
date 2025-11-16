using Random, ProgressMeter, Statistics

include("moves.jl")
include("helpers.jl")


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

                @showprogress 0.1 "runs" for run in 1:20
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
                    "lab3/$(data_name)/$(start_type)_$(local_search_func)_$(intra_function).png")
            end
        end
    end
end

const EdgeInfo = NamedTuple{(:orientation, :start_idx, :end_idx),Tuple{Symbol,Int,Int}}

struct SavedMove
    objective::Float64
    move_kind::Symbol
    old_edges::Vector{Tuple{Int,Int}}
    new_edges::Vector{Tuple{Int,Int}}
    orientations::Vector{Symbol}
end

edge_step_forward(idx, n) = idx == n ? 1 : idx + 1

function locate_edge(solution, edge)
    n = length(solution)
    a, b = edge

    idx_a = findfirst(isequal(a), solution)
    if !isnothing(idx_a)
        next_idx = edge_step_forward(idx_a, n)
        solution[next_idx] == b && return (:forward, idx_a, next_idx)
    end

    idx_b = findfirst(isequal(b), solution)
    if !isnothing(idx_b)
        next_idx = edge_step_forward(idx_b, n)
        solution[next_idx] == a && return (:reverse, idx_b, next_idx)
    end

    return (:missing, nothing, nothing)
end

function capture_edge_orientations(solution, old_edges)
    orientations = Symbol[]
    for edge in old_edges
        orientation, _, _ = locate_edge(solution, edge)
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

function apply_inter_move(solution, old_edges, new_edges)
    removed_node = find_repeated_node(old_edges)
    new_node = find_repeated_node(new_edges)

    if removed_node === nothing || new_node === nothing || in(new_node, solution)
        return nothing
    end

    idx = findfirst(isequal(removed_node), solution)
    isnothing(idx) && return nothing

    new_solution = copy(solution)
    new_solution[idx] = new_node
    return new_solution
end

function check_if_applicable(solution, move_kind, old_edges, new_edges, stored_orientations)
    length(old_edges) == length(stored_orientations) || return (:remove, nothing)

    edge_infos = EdgeInfo[]
    for edge in old_edges
        orientation, start_idx, end_idx = locate_edge(solution, edge)
        if orientation == :missing
            return (:remove, nothing)
        end
        push!(edge_infos, (orientation=orientation, start_idx=start_idx, end_idx=end_idx))
    end

    matches = true
    reversed = true
    for (idx, info) in enumerate(edge_infos)
        stored_orientation = stored_orientations[idx]
        matches &= info.orientation == stored_orientation
        reversed &= info.orientation != stored_orientation
    end

    if !(matches || reversed)
        return (:skip, nothing)
    end

    new_solution = move_kind == :intra ? apply_intra_move(solution, edge_infos) :
                   apply_inter_move(solution, old_edges, new_edges)

    new_solution === nothing && return (:remove, nothing)
    return (:apply, new_solution)
end


function steepest_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)
    improvement = true

    intra_move = [(i, j) for i in 1:n for j in i+2:n if abs(i - j) != n - 1]
    LM = SavedMove[]

    while improvement
        improvement = false

        # Evaluate all new moves and add improving moves to LM
        available_nodes = setdiff(1:size(distance_matrix, 1), solution)
        inter_move = [(i, j) for i in 1:n for j in available_nodes]

        moves = [(m, :intra) for m in intra_move]
        append!(moves, [(m, :inter) for m in inter_move])

        for ((a, b), move_kind) in moves
            if move_kind == :intra
                new_objective, old_edges, new_edges = intra_function(solution, objective, a, b, distance_matrix)
            else
                new_objective, old_edges, new_edges = inter_function(solution, objective, a, b, distance_matrix, costs)
            end
            if new_objective < objective
                orientations = capture_edge_orientations(solution, old_edges)
                orientations === nothing && continue
                push!(LM, SavedMove(Float64(new_objective), move_kind, old_edges, new_edges, orientations))
            end
        end

        sort!(LM, by=m -> m.objective, rev=true)

        i = length(LM)
        while i >= 1
            move = LM[i]
            status, new_solution = check_if_applicable(solution, move.move_kind, move.old_edges, move.new_edges, move.orientations)

            if status == :remove
                deleteat!(LM, i)
                i -= 1
            elseif status == :skip
                i -= 1
            elseif status == :apply
                refreshed_objective = calculate_cycle_length(new_solution, distance_matrix, costs)
                if refreshed_objective < objective
                    improvement = true
                    solution = new_solution
                    objective = refreshed_objective
                    deleteat!(LM, i)
                    break
                else
                    deleteat!(LM, i)
                    i -= 1
                end
            end
        end
    end

    return objective, solution
end