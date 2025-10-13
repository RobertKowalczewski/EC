using CSV, DataFrames, StatsBase, Plots

function euclidean(a, b)
    return round(sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2))
end

function create_distance_matrix(data_path)
    data = CSV.read(data_path, DataFrame; header=["x", "y", "w"])
    distance_matrix = zeros(nrow(data), nrow(data))

    for i in 1:nrow(data)
        for j in 1:nrow(data)
            if i != j
                distance_matrix[i, j] = euclidean(data[i, ["x", "y"]], data[j, ["x", "y"]]) + data[j, "w"]
            end
        end
    end

    return distance_matrix
end

# function find_best_solution(search_function, iterations, distance_matrix)
#     best_score = Inf
#     best_solution = nothing
#     current = 1

#     for i in 1:iterations
#         score, solution = search_function(distance_matrix, current)
#         current += 1
#         if score < best_score
#             best_score = score
#             best_solution = solution
#         end
#     end

#     return best_score, best_solution
# end

function find_best_solution_random(search_function, iterations, distance_matrix)
    scores = Float64[]
    solutions = Vector{Vector{Int}}(undef, iterations)
    for i in 1:iterations
        score, solution = search_function(distance_matrix)
        push!(scores, score)
        solutions[i] = solution
    end
    min_s = argmin(scores)
    max_s = argmax(scores)
    avg_s = mean(scores)
    println("Min: $(scores[min_s]) at run $(min_s)")
    println("Max: $(scores[max_s]) at run $(max_s)")
    println("Avg: $(avg_s)")
    return scores, solutions
end

function find_best_solution_greedy(search_function, distance_matrix)
    n = ceil(Int, size(distance_matrix)[1] / 2)
    scores = Float64[]
    solutions = Vector{Vector{Int}}(undef, n)
    for start in 1:n
        score, solution = search_function(distance_matrix, start)
        push!(scores, score)
        solutions[start] = solution
    end
    min_s = argmin(scores)
    max_s = argmax(scores)
    avg_s = mean(scores)
    println("Min: $(scores[min_s]) at start $(min_s)")
    println("Max: $(scores[max_s]) at start $(max_s)")
    println("Avg: $(avg_s)")
    return scores, solutions
end

function random_search(distance_matrix, start=nothing)
    n = ceil(Int, size(distance_matrix)[1] / 2)
    sampled_nodes = sample(1:size(distance_matrix)[1], n; replace=false)
    score = 0
    for i in 1:length(sampled_nodes)-1
        score += distance_matrix[sampled_nodes[i], sampled_nodes[i+1]]
    end

    score += distance_matrix[sampled_nodes[n], sampled_nodes[1]]
    return score, sampled_nodes
end

function nn_last_node(distance_matrix, start)
    #1. Choose a random node to start with
    n = size(distance_matrix)[1]
    cycle_length = ceil(Int, n / 2)

    solution = Vector{Int}(undef, cycle_length)
    visited = Set{Int}()

    score = 0
    current = start

    solution[1] = current
    push!(visited, current)

    for j in 2:cycle_length
        best_distance, best_node = Inf, nothing
        #2. Check all nodes for closest one including costs - iterate over row
        for i in 1:n
            if i ∉ visited
                distance = distance_matrix[current, i]
                if distance < best_distance
                    best_distance = distance
                    best_node = i
                end
            end
        end

        solution[j] = best_node
        push!(visited, best_node)
        score += best_distance
        current = best_node
    end

    score += distance_matrix[solution[end], solution[1]]

    return score, solution
end

# function nn_all(distance_matrix, start)
#     n = size(distance_matrix)[1]
#     cycle_length = ceil(Int, n / 2)

#     solution = Vector{Int}(undef, cycle_length)
#     visited = Set{Int}()

#     score = 0
#     current = start

#     solution[1] = current
#     push!(visited, current)

#     for node in 2:cycle_length
#         smallest_distance = Inf
#         smallest_distance_node = nothing
#         from = nothing # index after which to insert the new node
#         for i in 1:node
#             for j in 1:n
#                 distance = distance_matrix[solution[i], j]
#                 if distance < smallest_distance
#                     smallest_distance = distance
#                     smallest_distance_node = j
#                     from = i
#                 end
#             end
#         end
#     end

#     insert!()


# end


function nn_all(distance_matrix, start)
    n = size(distance_matrix)[1]
    cycle_length = ceil(Int, n / 2)
    solution = [start]
    visited = Set([start])

    # Add the second node: choose the nearest neighbor to start
    best_distance, best_node = Inf, nothing
    for i in 1:n
        if i ∉ visited
            d = distance_matrix[start, i]
            if d < best_distance
                best_distance = d
                best_node = i
            end
        end
    end
    push!(solution, best_node)
    push!(visited, best_node)

    # Insert remaining nodes at the best position
    while length(solution) < cycle_length
        # Find the closest unvisited node to any node in the current solution
        candidate, candidate_from = nothing, nothing
        min_dist = Inf
        for i in 1:n
            if i ∉ visited
                for j in solution
                    d = distance_matrix[j, i]
                    if d < min_dist
                        min_dist = d
                        candidate = i
                    end
                end
            end
        end

        # Find the best place to insert candidate in the cycle
        best_increase = Inf
        best_pos = 0
        for pos in 1:length(solution)
            prev = solution[pos]
            nxt = solution[mod1(pos+1, length(solution))]
            increase = distance_matrix[prev, candidate] + distance_matrix[candidate, nxt] - distance_matrix[prev, nxt]
            if increase < best_increase
                best_increase = increase
                best_pos = pos
            end
        end
        insert!(solution, best_pos+1, candidate)
        push!(visited, candidate)
    end

    # Calculate total score
    score = 0
    for i in 1:cycle_length
        score += distance_matrix[solution[i], solution[mod1(i+1, cycle_length)]]
    end

    return score, solution
end
k
# todo round mathematically euclidean
function greedy_cycle(distance_matrix, start=1)
    n = ceil(Int, size(distance_matrix)[1] / 2)
    solution = Array{Int}(undef, n)
    solution_set = Set{Int}()

    solution[1] = start
    push!(solution_set, start)

    #1. Choose closest node
    smallest_distance, closest_node = Inf, nothing
    for node in 1:size(distance_matrix)[1]
        if node != start
            d = distance_matrix[start, node]
            if d < smallest_distance
                smallest_distance = d
                closest_node = node
            end
        end
    end

    solution[2] = closest_node
    push!(solution_set, closest_node)

    for solution_node in 3:n
        best_cycle, best_node = Inf, nothing
        for node in 1:size(distance_matrix)[1]
            if node ∉ solution_set
                new_solution = Array{Int}(undef, solution_node)
                new_solution[1:solution_node-1] = solution[1:solution_node-1]
                new_solution[solution_node] = node

                cycle_length = calculate_cycle_length(new_solution, distance_matrix)
                if cycle_length < best_cycle
                    best_cycle = cycle_length
                    best_node = node
                end
            end
        end
        solution[solution_node] = best_node
        push!(solution_set, best_node)
    end
    return calculate_cycle_length(solution, distance_matrix), solution
end

# distance_matrix = create_distance_matrix("lab1/TSPA.csv")

distance_matrix = create_distance_matrix("lab1/TSPB.csv")


println("Running random search")
random_score, random_solution = find_best_solution_random(random_search, 200, distance_matrix)
#println("best score: ", random_score)
#println("solution: ", random_solution)

println("\nRunning nn_last search")
nn_score, nn_solution = find_best_solution_greedy(nn_last_node,distance_matrix)
#println("best score: ", nn_score)
#println("solution: ", nn_solution)

println("\nRunning nn_all search")
nn_all_score, nn_alll_solution = find_best_solution_greedy(nn_all,distance_matrix)
#println("best score: ", nn_all_score)
#println("solution: ", nn_alll_solution)

# best solution doesn't match for some reason? 
println("\nRunning greedy_cycle search")
greedy_score, greedy_solution = find_best_solution_greedy(greedy_cycle, distance_matrix)
#println("best score: ", greedy_score)





function report_best_solutions(method_results)
    best_scores = []
    best_solutions = []
    
    for (scores, solutions) in method_results
        min_idx = argmin(scores)
        push!(best_scores, scores[min_idx])
        push!(best_solutions, solutions[min_idx])
    end
    
    return best_scores, best_solutions
end

function print_zero_indexed_solutions(method_names, best_solutions)
    println("\n=== Best Solutions (0-indexed) ===")
    for (name, solution) in zip(method_names, best_solutions)
        # Convert from 1-indexed to 0-indexed
        zero_indexed = [node-1 for node in solution]
        println("$name:")
        println("[$(join(zero_indexed, ", "))]")
        println()

    end
    
    # Also print a consolidated list format for easy copying
    println("=== Solutions for Report ===")
    for (name, solution) in zip(method_names, best_solutions)
        zero_indexed = [node-1 for node in solution]
        println("$name: [$(join(zero_indexed, ", "))]")
    end
end

# After running all your search methods, collect results
method_results = [
    (random_score, random_solution),
    (nn_score, nn_solution),
    (nn_all_score, nn_alll_solution),
    (greedy_score, greedy_solution)
]

method_names = ["Random Search", "Nearest Neighbor", "NN All", "Greedy Cycle"]

# Find best solution for each method
best_scores, best_solutions = report_best_solutions(method_results)

# Print the best solutions with 0-indexed nodes
print_zero_indexed_solutions(method_names, best_solutions)

# You can also print a summary of the best scores
println("\n=== Best Scores ===")
for (name, score) in zip(method_names, best_scores)
    println("$name: $score")
end



# Modify the random_search to ensure it always includes node 1
function random_search_with_start(distance_matrix)
    n = ceil(Int, size(distance_matrix)[1] / 2)
    # Always include node 1
    sampled_nodes = [1]
    
    # Sample the remaining nodes
    remaining_nodes = setdiff(1:size(distance_matrix)[1], [1])
    needed_nodes = n - 1
    additional_nodes = sample(remaining_nodes, needed_nodes; replace=false)
    
    # Combine to form the complete solution
    append!(sampled_nodes, additional_nodes)
    
    score = 0
    for i in 1:length(sampled_nodes)-1
        score += distance_matrix[sampled_nodes[i], sampled_nodes[i+1]]
    end

    score += distance_matrix[sampled_nodes[n], sampled_nodes[1]]
    return score, sampled_nodes
end

# Run all methods again ensuring they start from node 1
println("\n=== Running all methods with node 1 as starting point ===")

println("\nRunning random search from node 1")
fixed_random_scores = Float64[]
fixed_random_solutions = Vector{Vector{Int}}()
for i in 1:100  # Run 100 times
    score, solution = random_search_with_start(distance_matrix)
    push!(fixed_random_scores, score)
    push!(fixed_random_solutions, solution)
end
min_idx = argmin(fixed_random_scores)
println("Min score: $(fixed_random_scores[min_idx])")

println("\nRunning nn_last from node 1")
nn_score_fixed, nn_solution_fixed = nn_last_node(distance_matrix, 1)
println("Score: $nn_score_fixed")

println("\nRunning nn_all from node 1")
nn_all_score_fixed, nn_all_solution_fixed = nn_all(distance_matrix, 1)
println("Score: $nn_all_score_fixed")

println("\nRunning greedy_cycle from node 1")
greedy_score_fixed, greedy_solution_fixed = greedy_cycle(distance_matrix, 1)
println("Score: $greedy_score_fixed")

# Collect fixed results that all start from node 1
fixed_method_results = [
    (fixed_random_scores[min_idx], fixed_random_solutions[min_idx]),
    (nn_score_fixed, nn_solution_fixed),
    (nn_all_score_fixed, nn_all_solution_fixed),
    (greedy_score_fixed, greedy_solution_fixed)
]

# Pretty print the solutions for the report
println("\n=== Best Solutions (0-indexed) starting from node 0 ===")
for (name, (score, solution)) in zip(method_names, fixed_method_results)
    # Convert from 1-indexed to 0-indexed
    zero_indexed = [node-1 for node in solution]
    println("$name: [$(join(zero_indexed, ", "))]")
    println("Score: $score")
    println()
end

# Print a consolidated list format for easy copying
println("=== Solutions for Report (all starting from node 0) ===")
for (name, (score, solution)) in zip(method_names, fixed_method_results)
    zero_indexed = [node-1 for node in solution]
    println("$name: [$(join(zero_indexed, ", "))]")
end