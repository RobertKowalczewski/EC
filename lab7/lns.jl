include("local-search.jl")
include("greedy-cycle-2-regret.jl")
include("destroy.jl")
include("helpers.jl")


function lns(distance_matrix, costs, run_time, use_local_search)
    x = random_start(distance_matrix)

    x_objective, x = steepest_local_search(x, distance_matrix, costs)
    start_time = time()
    i = 0

    while (time() - start_time) <= run_time
        i += 1
        y = copy(x)
        destroy!(y)
        if use_local_search
            y = greedy_cycle_2_regret(distance_matrix, costs, [0.5, 0.5], y)
            y_objective, y = steepest_local_search(y, distance_matrix, costs)
        else
            y = greedy_cycle_2_regret(distance_matrix, costs, [0.5, 0.5], y)
            y_objective = calculate_cycle_length(y, distance_matrix, costs)
        end
        if y_objective < x_objective
            x = y
            x_objective = y_objective
        end
    end
    return x_objective, x, i
end