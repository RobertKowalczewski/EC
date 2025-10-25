using Random, ProgressMeter

include("moves.jl")
include("helpers.jl")
include("nn-greedy-2-regret.jl")


function test_algorithms(distance_matrix, costs)
    @showprogress 0.5 for start_type in [random_start, nn_greedy_2_regret
        ] # TODO
        @showprogress 0.5 for local_search_func in [greedy_local_search] # TODO
            @showprogress 0.5 for intra_function in [intra_two_nodes_exchange]#, intra_two_edges_exchange] # TODO
                
                objectives = []
                solutions = []

                @showprogress 0.5 "runs" for run in 1:200
                    starting_solution = start_type(distance_matrix,costs, [0.5, 0.5], run)
                    objective, solution = local_search_func(starting_solution, distance_matrix, costs, intra_function, inter_two_nodes_exchange)
                    push!(objectives, objective)
                    push!(solutions, solution)
                end

                min_s = argmin(objectives)
                max_s = argmax(objectives)
                avg_s = mean(objectives)
                println(start_type," ", local_search_func," ", intra_function)
                println("Min: $(objectives[min_s]) at start $(min_s)")
                println("Max: $(objectives[max_s]) at start $(max_s)")
                println("Avg: $(avg_s)")
            end
        end
    end
end

function greedy_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)
    improvement = true

    intra_move = [(i, j) for i in 1:n for j in i+1:n]

    while improvement
        improvement = false

        available_nodes = setdiff(1:size(distance_matrix, 1), solution)
        inter_move = [(i, j) for i in 1:n for j in available_nodes]

        moves = [(m, :intra) for m in intra_move]
        append!(moves, [(m, :inter) for m in inter_move])
        shuffle!(moves)


        for ((a,b), move_kind) in moves
            if move_kind == :intra
                new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
            else
                new_objective, new_solution = Inf, nothing#inter_function(solution, objective, a, b, distance_matrix, costs)
            end
            if new_objective < objective

                solution = new_solution
                objective = new_objective
                improvement = true
                break
            end
        end
    end

    return objective, solution
end

function steepest_local_search(starting_solution, distance_matrix, costs, intra_function, inter_function)
    n = length(starting_solution)
    solution = starting_solution
    objective = calculate_cycle_length(starting_solution, distance_matrix, costs)
    improvement = true

    intra_move = [(i, j) for i in 1:n for j in i+1:n]

    while improvement
        improvement = false

        available_nodes = setdiff(1:size(distance_matrix, 1), solution)
        inter_move = [(i, j) for i in 1:n for j in available_nodes]

        moves = [(m, :intra) for m in intra_move]
        append!(moves, [(m, :inter) for m in inter_move])
        shuffle!(moves)


        for ((a,b), move_kind) in moves
            if move_kind == :intra
                new_objective, new_solution = intra_function(solution, objective, a, b, distance_matrix)
            else
                new_objective, new_solution = Inf, nothing#inter_function(solution, objective, a, b, distance_matrix, costs)
            end
            if new_objective < objective

                solution = new_solution
                objective = new_objective
                improvement = true
                break
            end
        end
    end

    return objective, solution
end
