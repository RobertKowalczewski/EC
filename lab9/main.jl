using Random
using StatsBase
using Statistics

include("helpers.jl")

# --- Move utilities ---------------------------------------------------------

@inline function intra_two_edges_delta(solution, a, b, distance_matrix, n)
	a2 = a == n ? 1 : a + 1
	b2 = b == n ? 1 : b + 1

	old_edges = distance_matrix[solution[a], solution[a2]] + distance_matrix[solution[b], solution[b2]]
	new_edges = distance_matrix[solution[a], solution[b]] + distance_matrix[solution[a2], solution[b2]]
	return new_edges - old_edges
end

@inline function inter_two_nodes_delta(solution, a, new_node, distance_matrix, costs, n)
	prev_a = a == 1 ? n : a - 1
	next_a = a == n ? 1 : a + 1
	current_node = solution[a]

	old_edges = distance_matrix[solution[prev_a], current_node] + distance_matrix[current_node, solution[next_a]]
	new_edges = distance_matrix[solution[prev_a], new_node] + distance_matrix[new_node, solution[next_a]]
	return (new_edges - old_edges) + (costs[new_node] - costs[current_node])
end

@inline function apply_intra_two_edges!(solution, a, b)
	left = a + 1
	right = b
	while left < right
		solution[left], solution[right] = solution[right], solution[left]
		left += 1
		right -= 1
	end
end

@inline function apply_inter_two_nodes!(solution, a, new_node)
	old_node = solution[a]
	solution[a] = new_node
	return old_node
end

# --- Local search -----------------------------------------------------------

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

# --- Construction heuristics ------------------------------------------------

function greedy_cycle_2_regret(distance_matrix, costs, weights, incomplete_solution)
	n = size(distance_matrix, 1)
	cycle_length = ceil(Int, n / 2)
	solution = copy(incomplete_solution)
	visited = Set(solution)

	while length(solution) < cycle_length
		best_score = -Inf
		best_node = nothing
		best_insertion_pos = nothing

		for i in 1:n
			if i ∉ visited
				best_cost = Inf
				second_best_cost = Inf
				best_pos = nothing

				for pos in eachindex(solution)
					prev = solution[pos]
					nxt = solution[mod1(pos + 1, length(solution))]
					insertion_cost = distance_matrix[prev, i] + costs[i] + distance_matrix[i, nxt] - distance_matrix[prev, nxt]

					if insertion_cost < best_cost
						second_best_cost = best_cost
						best_cost = insertion_cost
						best_pos = pos + 1
					elseif insertion_cost < second_best_cost
						second_best_cost = insertion_cost
					end
				end

				regret = second_best_cost - best_cost
				score = weights[1] * regret - weights[2] * best_cost

				if score > best_score
					best_score = score
					best_node = i
					best_insertion_pos = best_pos
				end
			end
		end

		insert!(solution, best_insertion_pos, best_node)
		push!(visited, best_node)
	end

	return solution
end

function destroy!(x, p=0.3)
	remove_indexes = sort(sample(1:length(x), floor(Int, length(x) * p), replace=false))
	deleteat!(x, remove_indexes)
end

# --- Evolutionary operators -------------------------------------------------

solution_key(sol) = join(sol, ",")

function extract_common_subpaths(p1, p2, total_nodes)
	# Identify common edges (undirected) and return vector of subpaths; also capture common nodes without shared edges
	edges1 = Set{Tuple{Int,Int}}()
	edges2 = Set{Tuple{Int,Int}}()
	n1 = length(p1)
	n2 = length(p2)
	for i in 1:n1
		a, b = p1[i], p1[mod1(i + 1, n1)]
		push!(edges1, a < b ? (a, b) : (b, a))
	end
	for i in 1:n2
		a, b = p2[i], p2[mod1(i + 1, n2)]
		push!(edges2, a < b ? (a, b) : (b, a))
	end

	common_edges = intersect(edges1, edges2)
	common_nodes = intersect(Set(p1), Set(p2))

	# If no common structure at all, return empty and let caller randomize
	if isempty(common_edges) && isempty(common_nodes)
		return Vector{Vector{Int}}()
	end

	adjacency = [Int[] for _ in 1:total_nodes]
	for (a, b) in common_edges
		push!(adjacency[a], b)
		push!(adjacency[b], a)
	end

	visited = falses(total_nodes)
	subpaths = Vector{Vector{Int}}()

	for node in 1:total_nodes
		if !visited[node] && !isempty(adjacency[node])
			if length(adjacency[node]) == 1
				# Open path start
				current = node
				path = Int[]
				prev = 0
				while true
					push!(path, current)
					visited[current] = true
					next_nodes = [n for n in adjacency[current] if n != prev]
					if isempty(next_nodes)
						break
					end
					prev, current = current, next_nodes[1]
					if visited[current]
						push!(path, current)
						break
					end
				end
				push!(subpaths, unique(path))
			end
		end
	end

	# Handle potential cycles composed only of degree-2 nodes
	for node in 1:total_nodes
		if !visited[node] && !isempty(adjacency[node])
			path = Int[]
			current = node
			prev = 0
			while true
				push!(path, current)
				visited[current] = true
				next_nodes = [n for n in adjacency[current] if n != prev]
				if isempty(next_nodes)
					break
				end
				prev, current = current, next_nodes[1]
				if visited[current]
					break
				end
			end
			push!(subpaths, unique(path))
		end
	end

	# Add singleton subpaths for common nodes that are not part of any common edge
	for node in common_nodes
		if !visited[node]
			push!(subpaths, [node])
			visited[node] = true
		end
	end

	return subpaths
end

function recombine_operator1(p1, p2, total_nodes)
	target_len = ceil(Int, total_nodes / 2)
	subpaths = extract_common_subpaths(p1, p2, total_nodes)
	used = Set{Int}()
	for path in subpaths
		foreach(node -> push!(used, node), path)
	end

	# Add random singleton subpaths to reach target size
	remaining_nodes = [n for n in 1:total_nodes if n ∉ used]
	shuffle!(remaining_nodes)
	for node in remaining_nodes
		if length(used) >= target_len
			break
		end
		push!(subpaths, [node])
		push!(used, node)
	end

	# Randomly connect subpaths
	shuffle!(subpaths)
	if isempty(subpaths)
		return randperm(total_nodes)[1:target_len]
	end

	route = copy(subpaths[1])
	for idx in 2:length(subpaths)
		path = subpaths[idx]
		if rand(Bool)
			reverse!(path)
		end
		if rand(Bool)
			route = vcat(route, path)
		else
			route = vcat(path, route)
		end
	end

	return route[1:target_len]
end

function recombine_operator2(p1, p2, distance_matrix, costs)
	target_len = ceil(Int, size(distance_matrix, 1) / 2)
	base, other = rand(Bool) ? (p1, p2) : (p2, p1)
	other_set = Set(other)

	shared = [node for node in base if node in other_set]
	shared = unique(shared)
	if isempty(shared)
		shared = [rand(1:length(other_set))]
	end

	offspring = greedy_cycle_2_regret(distance_matrix, costs, [0.5, 0.5], shared)

	return offspring[1:target_len]
end

# --- Metaheuristics ---------------------------------------------------------

function build_initial_solution(distance_matrix, costs)
	start = random_start(distance_matrix)
	greedy = greedy_cycle_2_regret(distance_matrix, costs, [0.5, 0.5], start)
	_, best = steepest_local_search(greedy, distance_matrix, costs)
	return best
end

function build_population(pop_size, distance_matrix, costs)
	population = Vector{Vector{Int}}()
	scores = Float64[]
	seen = Set{String}()

	while length(population) < pop_size
		sol = build_initial_solution(distance_matrix, costs)
		key = solution_key(sol)
		if key in seen
			continue
		end
		push!(population, sol)
		push!(scores, calculate_cycle_length(sol, distance_matrix, costs))
		push!(seen, key)
	end

	return population, scores, seen
end

function steady_state_hybrid(distance_matrix, costs; pop_size=20, time_limit=10.0, operator=:op1, apply_ls_offspring=true)
	start_time = time()
	total_nodes = size(distance_matrix, 1)
	population, scores, seen = build_population(pop_size, distance_matrix, costs)

	best_score, best_idx = findmin(scores)
	best_solution = copy(population[best_idx])

	while time() - start_time < time_limit
		p_indices = sample(1:pop_size, 2; replace=false)
		p1 = population[p_indices[1]]
		p2 = population[p_indices[2]]

		child = operator == :op1 ? recombine_operator1(p1, p2, total_nodes) : recombine_operator2(p1, p2, distance_matrix, costs)
		child_score = calculate_cycle_length(child, distance_matrix, costs)

		if apply_ls_offspring
			child_score, child = steepest_local_search(child, distance_matrix, costs)
		end

		key = solution_key(child)
		if key in seen
			continue
		end

		worst_score, worst_idx = findmax(scores)
		if child_score < worst_score
			# Replace worst individual
			delete!(seen, solution_key(population[worst_idx]))
			population[worst_idx] = child
			scores[worst_idx] = child_score
			push!(seen, key)

			if child_score < best_score
				best_score = child_score
				best_solution = copy(child)
			end
		end
	end

	return best_score, best_solution
end

# --- Driver -----------------------------------------------------------------

function run_instance(data_path; time_limit=10.0)
	distance_matrix, costs = prepare_data(data_path)

	function summarize_runs(operator; apply_ls_offspring=true, label="")
		scores = Float64[]
		solutions = Vector{Vector{Int}}()
		best_score = Inf
		best_solution = Int[]
		for _ in 1:20
			score, sol = steady_state_hybrid(distance_matrix, costs; pop_size=20, time_limit=time_limit, operator=operator, apply_ls_offspring=apply_ls_offspring)
			push!(scores, score)
			push!(solutions, sol)
			if score < best_score
				best_score = score
				best_solution = sol
			end
		end
		return (label=label, avg=mean(scores), min=minimum(scores), max=maximum(scores), scores=scores, solutions=solutions, best_score=best_score, best_solution=best_solution)
	end

	stats = [
		summarize_runs(:op1; apply_ls_offspring=true, label="EA op1"),
		summarize_runs(:op2; apply_ls_offspring=true, label="EA op2 LS"),
		summarize_runs(:op2; apply_ls_offspring=false, label="EA op2 no LS")
	]

	for s in stats
		println(string(s.label, ": ", s.avg, " (", s.min, ", ", s.max, ")"))
	end

	for s in stats
		title = s.label * " $(data_path)"
		save_name = "lab9/plots/" * replace(s.label, ' ' => "_") * "_$(data_path)_best.png"
		plot_best_solution(s.scores, s.solutions, data_path, title, save_name)
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	data_path = length(ARGS) >= 1 ? ARGS[1] : "TSPB.csv"
	if data_path == "TSPA.csv"
		time_limit = 1.59246505045
	else
		time_limit = 1.9890074519000003
	end
	Random.seed!(1234)
	run_instance(data_path, time_limit=time_limit)
end
