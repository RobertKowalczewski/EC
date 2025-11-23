using Random, ProgressMeter, Statistics

include("helpers.jl")
include("moves.jl")

function test_algorithms(distance_matrix, costs, data_path, data_name)
    for start_type in [random_start]
        for local_search_func in [steepest_local_search]
            for intra_function in [intra_two_nodes_exchange]

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

function perturb(x)
end

function iterated_local_search(distance_matrix, costs, run_time)
    x = random_start(distance_matrix)
    objective, x = steepest_local_search(x, distance_matrix, costs, intra_two_edges_exchange, inter_two_nodes_exchange)
    t = time()
    while (time() - t) <= run_time
        y = perturb(x)
        objective_y, y = steepest_local_search(y, distance_matrix, costs, intra_two_edges_exchange, inter_two_nodes_exchange)
        if objective_y < objective
            objective, x = objective_y, y
        end
    end
    return objective, x
end

function steepest_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)
    improvement = true

    if intra_function == intra_two_nodes_exchange
        intra_move = [(i, j) for i in 1:n for j in i+1:n]
    elseif intra_function == intra_two_edges_exchange
        intra_move = [(i, j) for i in 1:n for j in i+2:n if abs(i - j) != n - 1]
    else
        println("invalide intra function")
        return
    end

    while improvement
        improvement = false

        available_nodes = setdiff(1:size(distance_matrix, 1), solution)
        inter_move = [(i, j) for i in 1:n for j in available_nodes]

        moves = [(m, :intra) for m in intra_move]
        append!(moves, [(m, :inter) for m in inter_move])

        best_objective, best_solution = Inf, nothing
        for ((a, b), move_kind) in moves
            if move_kind == :intra
                new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
            else
                new_objective, new_solution = inter_function(solution, objective, a, b, distance_matrix, costs)
            end
            if new_objective < objective && new_objective < best_objective
                best_solution = new_solution
                best_objective = new_objective
                improvement = true
            end
        end

        if improvement
            solution, objective = best_solution, best_objective
        end
    end

    return objective, solution
end