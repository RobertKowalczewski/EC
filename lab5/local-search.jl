using Random, ProgressMeter, Statistics

include("moves.jl")
include("helpers.jl")
include("nn-greedy-2-regret.jl")
include("greedy-cycle-2-regret.jl")


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
                    "lab3/$(data_name)/$(start_type)_$(local_search_func)_$(intra_function).png")
            end
        end
    end
end

function check_if_applicable(solution, move_kind, old_edges, new_edges)
    if move_kind == :inter
        # Removed edges (defining the saved move) no longer exist in the current solution (at least one of them)
        old_edge_removed = false
        for (a, b) in old_edges
            check_edge_exists(solution, a, b)
        end

    else

    end
end

function steepest_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)
    improvement = true

    intra_move = [(i, j) for i in 1:n for j in i+2:n if abs(i - j) != n - 1]
    LM = []

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
                append!(LM, [new_objective, move_kind, old_edges, new_edges])
            end
        end

        sort!(LM, by=x -> x[1], rev=true)

        # For moves m from LM starting from the best until a applicable move is found
        for (new_objective, move_kind, old_edges, new_edges) in LM
            # Check if m is applicable and if not remove it from LM
            applicable, new_solution = check_if_applicable(solution, move_kind, old_edges, new_edges)
            if applicable
                improvement = true
                solution = new_solution
                objective = new_objective # ?
                break
            end
        end




    end

    return objective, solution
end