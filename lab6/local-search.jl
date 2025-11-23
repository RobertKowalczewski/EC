using Random, ProgressMeter, Statistics

include("helpers.jl")
include("moves.jl")

function test_algorithms(distance_matrix, costs, data_path, data_name)
    objectives = []
    solutions = []
    starts = []
    start_objectives = []
    times = []

    fmt(x) = round(x; digits=5)

    time = 10

    @showprogress 0.1 "runs" for run in 1:20
        starting_solution = random_start(distance_matrix)
        t = @elapsed begin
            objective, solution = iterated_local_search(distance_matrix, costs, time)
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

"""Apply random swaps and node replacements to diversify the route."""
function perturb(solution, distance_matrix, perturb_size)
    perturbed = copy(solution)
    n = length(solution)
    total_nodes = size(distance_matrix, 1)
    available_nodes = Set(setdiff(1:total_nodes, solution))
    steps = max(1, Int(ceil(perturb_size)))

    for _ in 1:steps
        perform_swap = isempty(available_nodes) ? true : rand() < 0.5

        if perform_swap && n > 1
            i = rand(1:n)
            j = rand(1:n)
            while j == i
                j = rand(1:n)
            end
            perturbed[i], perturbed[j] = perturbed[j], perturbed[i]
        else
            idx = rand(1:n)
            new_node = rand(collect(available_nodes))
            old_node = perturbed[idx]
            delete!(available_nodes, new_node)
            push!(available_nodes, old_node)
            perturbed[idx] = new_node
        end
    end

    return perturbed
end

function iterated_local_search(distance_matrix, costs, run_time; perturb_size=2)
    x = random_start(distance_matrix)
    objective, x = steepest_local_search(x, distance_matrix, costs, intra_two_edges_exchange, inter_two_nodes_exchange)
    start_time = time()

    while (time() - start_time) <= run_time
        y = perturb(x, distance_matrix, perturb_size)
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

    intra_move = [(i, j) for i in 1:n for j in i+2:n if abs(i - j) != n - 1]

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