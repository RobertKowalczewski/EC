using Random, ProgressMeter, Statistics
using Base.Threads

include("helpers.jl")
include("moves.jl")

const NUM_RUNS = 20

fmt(x) = round(x; digits=5)

"""Run independent starts in parallel and capture their statistics."""
function run_parallel_runs(num_runs, distance_matrix, costs, run_solver; record_iterations::Bool=false)
    objectives = Vector{Float64}(undef, num_runs)
    solutions = Vector{Vector{Int}}(undef, num_runs)
    starts = Vector{Vector{Int}}(undef, num_runs)
    start_objectives = Vector{Float64}(undef, num_runs)
    times = Vector{Float64}(undef, num_runs)
    iter_counts = record_iterations ? Vector{Int}(undef, num_runs) : nothing

    progress = Progress(num_runs; desc="runs", dt=0.1)
    progress_lock = SpinLock()

    @threads for run in 1:num_runs
        rng = Random.TaskLocalRNG()
        start_solution = random_start(distance_matrix; rng=rng)
        starts[run] = copy(start_solution)
        start_objectives[run] = calculate_cycle_length(start_solution, distance_matrix, costs)

        elapsed_time = @elapsed begin
            result = run_solver(start_solution, rng)
            if record_iterations
                objective, solution, iterations = result
                iter_counts[run] = iterations
            else
                objective, solution = result
            end
            objectives[run] = Float64(objective)
            solutions[run] = solution
        end

        times[run] = elapsed_time

        lock(progress_lock) do
            ProgressMeter.next!(progress)
        end
    end

    return (; objectives, solutions, starts, start_objectives, times, iterations=iter_counts)
end

function summarize_runs(run_data)
    objectives = run_data.objectives
    start_objectives = run_data.start_objectives
    times = run_data.times

    min_idx = argmin(objectives)
    max_idx = argmax(objectives)

    return (
        avg_objective=mean(objectives),
        min_objective=objectives[min_idx],
        max_objective=objectives[max_idx],
        avg_start=mean(start_objectives),
        min_start=minimum(start_objectives),
        max_start=maximum(start_objectives),
        avg_time=mean(times),
        min_time=minimum(times),
        max_time=maximum(times),
        best_index=min_idx,
    )
end

function test_algorithms(distance_matrix, costs, data_path, data_name)
    println("MSLS:")
    msls_runs = run_parallel_runs(
        NUM_RUNS,
        distance_matrix,
        costs,
        (start_solution, rng) -> begin
            multiple_start_local_search(distance_matrix, costs; rng=rng)
        end,
    )
    msls_stats = summarize_runs(msls_runs)
    report_algorithm(msls_stats)
    plot_best_solution(msls_runs.objectives, msls_runs.solutions,
        data_path,
        "$(data_name)_MSLS",
        "$(data_name)_multiple_start_LS.png")

    println("solution from first node:")
    println([x - 1 for x in msls_runs.solutions[1]])

    run_time = max(msls_stats.avg_time, eps())

    println("ILS:")
    ils_runs = run_parallel_runs(
        NUM_RUNS,
        distance_matrix,
        costs,
        (start_solution, rng) -> begin
            iterated_local_search(distance_matrix, costs, run_time;
                rng=rng,
                initial_solution=copy(start_solution))
        end;
        record_iterations=true,
    )
    ils_stats = summarize_runs(ils_runs)
    report_algorithm(ils_stats)
    println("number of basic iterations of LS: $(ils_runs.iterations)")
    plot_best_solution(ils_runs.objectives, ils_runs.solutions,
        data_path,
        "$(data_name)_ILS",
        "$(data_name)_iterated_local_search.png")

    println("solution from first node:")
    println([x - 1 for x in ils_runs.solutions[1]])
end

function report_algorithm(stats)
    println("After Local Search | Starting Solution")
    println("Obj: $(fmt(stats.avg_objective)) ($(fmt(stats.min_objective)) - $(fmt(stats.max_objective))) | Obj: $(fmt(stats.avg_start)) ($(fmt(stats.min_start)) - $(fmt(stats.max_start)))")
    println("Time[s]: $(fmt(stats.avg_time)) ($(fmt(stats.min_time)) - $(fmt(stats.max_time)))")
end

"""Apply random swaps and node replacements to diversify the route."""
function perturb(solution, distance_matrix, perturb_size; rng=Random.default_rng())
    perturbed = copy(solution)
    n = length(solution)
    total_nodes = size(distance_matrix, 1)
    steps = max(1, Int(ceil(perturb_size)))

    available_nodes = Vector{Int}()
    in_solution = falses(total_nodes)
    @inbounds for node in perturbed
        in_solution[node] = true
    end
    for node in 1:total_nodes
        if !in_solution[node]
            push!(available_nodes, node)
        end
    end

    for _ in 1:steps
        perform_swap = isempty(available_nodes) ? true : rand(rng) < 0.5

        if perform_swap && n > 1
            i = rand(rng, 1:n)
            j = rand(rng, 1:n)
            while j == i
                j = rand(rng, 1:n)
            end
            perturbed[i], perturbed[j] = perturbed[j], perturbed[i]
        else
            idx = rand(rng, 1:n)
            avail_idx = rand(rng, 1:length(available_nodes))
            new_node = available_nodes[avail_idx]
            old_node = perturbed[idx]
            perturbed[idx] = new_node
            available_nodes[avail_idx] = old_node
        end
    end

    return perturbed
end

function iterated_local_search(distance_matrix, costs, run_time; perturb_size=2, rng=Random.default_rng(), initial_solution=nothing)
    x = isnothing(initial_solution) ? random_start(distance_matrix; rng=rng) : copy(initial_solution)
    objective, x = steepest_local_search(x, distance_matrix, costs)
    start_time = time()
    i = 0

    while (time() - start_time) <= run_time
        i += 1
        y = perturb(x, distance_matrix, perturb_size; rng=rng)
        objective_y, y = steepest_local_search(y, distance_matrix, costs)
        if objective_y < objective
            objective, x = objective_y, y
        end
    end

    return objective, x, i
end

function multiple_start_local_search(distance_matrix, costs; rng=Random.default_rng())
    objective = Inf
    solution = Vector{Int}()
    for _ in 1:200
        start_solution = random_start(distance_matrix; rng=rng)
        new_objective, new_solution = steepest_local_search(start_solution, distance_matrix, costs)

        if new_objective < objective
            objective = new_objective
            solution = new_solution
        end
    end

    return objective, solution
end

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