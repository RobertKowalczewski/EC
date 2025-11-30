using ProgressMeter

include("lns.jl")

function test_lns(distance_matrix, costs, data_name, run_time, data_path)
    n_runs = 20

    objectives = []
    iterations = []
    solutions = []

    @showprogress 0.1 "runs" for i in 1:n_runs
        x_objective, x, i = lns(distance_matrix, costs,run_time,  true)
        append!(objectives, x_objective)
        append!(iterations, i)
        push!(solutions, x)
    end

    println("LNS with LS")
    min_obj = argmin(objectives)
    max_obj = argmax(objectives)
    avg_obj = mean(objectives)

    min_i = argmin(iterations)
    max_i = argmax(iterations)
    avg_i = mean(iterations)

    println("Obj: $(fmt(avg_obj)) ($(fmt(objectives[min_obj])) - $(fmt(objectives[max_obj])))")
    println("Iterations: $(avg_i) ($(iterations[min_i]) - $(iterations[max_i]))")

    plot_best_solution(
        objectives, 
        solutions,
        data_path,
        "large neighborhood search with LS - $(data_name)",
        "lab7/$(data_name)/lns_ls.png"
    )

    objectives = []
    iterations = []
    solutions = []

    @showprogress 0.1 "runs" for i in 1:n_runs
        x_objective, x, i = lns(distance_matrix, costs,run_time, false)
        append!(objectives, x_objective)
        append!(iterations, i)
        push!(solutions, x)
    end

    println("LNS without LS")
    min_obj = argmin(objectives)
    max_obj = argmax(objectives)
    avg_obj = mean(objectives)

    min_i = argmin(iterations)
    max_i = argmax(iterations)
    avg_i = mean(iterations)

    println("Obj: $(fmt(avg_obj)) ($(fmt(objectives[min_obj])) - $(fmt(objectives[max_obj])))")
    println("Iterations: $(avg_i) ($(iterations[min_i]) - $(iterations[max_i]))")

    plot_best_solution(
        objectives, 
        solutions,
        data_path,
        "large neighborhood search without LS - $(data_name)",
        "lab7/$(data_name)/lns.png"
    )
end

println("TSPA:")
run_time = 0.61072
data_path = "lab1/TSPA.csv"
distance_matrix, costs = prepare_data(data_path)

test_lns(distance_matrix, costs, "TSPA", run_time, data_path)

println("TSPB:")
run_time = 0.61414
data_path = "lab1/TSPB.csv"
distance_matrix, costs = prepare_data(data_path)

test_lns(distance_matrix, costs, "TSPB", run_time, data_path)