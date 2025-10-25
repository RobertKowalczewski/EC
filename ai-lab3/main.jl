using Random
using Statistics
using Printf
using Base.Filesystem: basename

include("helpers.jl")
include("local-search.jl")

const DEFAULT_RUNS = 200
const DEFAULT_SEED = 1_234_567

function resolve_dataset_path(path::AbstractString)
	return isabspath(path) ? path : normpath(joinpath(@__DIR__, path))
end

function default_dataset_paths()
	return [normpath(joinpath(@__DIR__, "..", "lab1", "TSPA.csv")),
			normpath(joinpath(@__DIR__, "..", "lab1", "TSPB.csv"))]
end

function greedy_start_sequence(instance::TSPInstance, runs::Int, rng::AbstractRNG)
	n = length(instance)
	order = collect(1:n)
	shuffle!(rng, order)
	sequence = Vector{Int}(undef, runs)
	for idx in 1:runs
		pos = mod1(idx, n)
		sequence[idx] = order[pos]
		if pos == n && idx < runs
			shuffle!(rng, order)
		end
	end
	return sequence
end

function run_configuration(instance::TSPInstance, local_type::Symbol, intra_type::Symbol, start_type::Symbol; runs::Int=DEFAULT_RUNS, seed::Int=DEFAULT_SEED, verify::Bool=false)
	rng = MersenneTwister(seed)
	objectives = Vector{Int}(undef, runs)
	initial_objectives = Vector{Int}(undef, runs)
	times = Vector{Float64}(undef, runs)
	iterations = Vector{Int}(undef, runs)
	evaluated_moves = Vector{Int}(undef, runs)
	applied_moves = Vector{Int}(undef, runs)
	start_nodes = start_type == :greedy ? greedy_start_sequence(instance, runs, rng) : nothing

	for run in 1:runs
		initial_solution = start_type == :random ? random_initial_solution(instance, rng) : two_regret_greedy_initial_solution(instance, start_nodes[run])
		initial_objectives[run] = objective(initial_solution)
		stats = SearchStats()
		timed = @timed begin
			if local_type == :steepest
				steepest_local_search(initial_solution, instance, intra_type; stats=stats)
			else
				greedy_local_search(initial_solution, instance, intra_type; rng=rng, stats=stats)
			end
		end
		final_solution, stats = timed.value
		if verify
			expected = recompute_objective(final_solution, instance)
			@assert expected == final_solution.objective "Objective cache mismatch: expected $(expected), got $(final_solution.objective)"
		end
		objectives[run] = final_solution.objective
		times[run] = timed.time
		iterations[run] = stats.iterations
		evaluated_moves[run] = stats.evaluated_moves
		applied_moves[run] = stats.applied_moves
	end

	improvements = initial_objectives .- objectives
	return (; objectives, initial_objectives, times, iterations, evaluated_moves, applied_moves, improvements)
end

function summarize_configuration(data)
	obj = Float64.(data.objectives)
	init = Float64.(data.initial_objectives)
	improvements = Float64.(data.improvements)
	iterations = Float64.(data.iterations)
	evaluated = Float64.(data.evaluated_moves)
	applied = Float64.(data.applied_moves)
	return (; best = minimum(data.objectives),
			worst = maximum(data.objectives),
			mean = mean(obj),
			std = std(obj),
			mean_initial = mean(init),
			mean_improvement = mean(improvements),
			best_improvement = maximum(data.improvements),
			mean_iterations = mean(iterations),
			mean_evaluated = mean(evaluated),
			mean_applied = mean(applied),
			mean_time = mean(data.times))
end

function print_configuration_summary(dataset_label::AbstractString, local_type::Symbol, intra_type::Symbol, start_type::Symbol, summary)
	@printf("    %-8s %-8s %-8s | best=%8d mean=%10.2f std=%10.2f init_mean=%10.2f mean_impr=%9.2f best_impr=%6d iter=%7.2f eval=%10.2f appl=%10.2f time=%7.4f\n",
		string(local_type), string(intra_type), string(start_type), summary.best, summary.mean, summary.std,
		summary.mean_initial, summary.mean_improvement, summary.best_improvement, summary.mean_iterations,
		summary.mean_evaluated, summary.mean_applied, summary.mean_time)
end

function run_experiments(dataset_paths::Vector{String}; runs::Int=DEFAULT_RUNS, seed::Int=DEFAULT_SEED, verify::Bool=false)
	local_types = [:steepest, :greedy]
	intra_types = [:node, :edge]
	start_types = [:random, :greedy]
	config_index = 0

	for dataset in dataset_paths
		resolved = resolve_dataset_path(dataset)
		instance = load_instance(resolved)
		dataset_label = basename(resolved)
		println("Dataset: $(dataset_label) | nodes=$(length(instance)) | required=$(instance.required)")

		for local_type in local_types
			for intra_type in intra_types
				for start_type in start_types
					config_index += 1
					config_seed = seed + config_index
					data = run_configuration(instance, local_type, intra_type, start_type; runs=runs, seed=config_seed, verify=verify)
					summary = summarize_configuration(data)
					print_configuration_summary(dataset_label, local_type, intra_type, start_type, summary)
				end
			end
		end

		println()
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	dataset_args = isempty(ARGS) ? default_dataset_paths() : ARGS
	run_experiments(collect(dataset_args))
end
