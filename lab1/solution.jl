using CSV, DataFrames, StatsBase

function euclidean(a, b)
    return sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
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

function find_best_solution(search_function, iterations, distance_matrix)
    best_score = Inf
    best_solution = nothing
    current = 1

    for i in 1:iterations
        score, solution = search_function(distance_matrix, current)
        current += 1
        if score < best_score
            best_score = score
            best_solution = solution
        end
    end

    return best_score, best_solution
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

function nn_all(distance_matrix, start)
    n = size(distance_matrix)[1]
    cycle_length = ceil(Int, n / 2)

    solution = Vector{Int}(undef, cycle_length)
    visited = Set{Int}()

    score = 0
    current = start

    solution[1] = current
    push!(visited, current)

    for node in 2:cycle_length
        smallest_distance = Inf
        smallest_distance_node = nothing
        from = nothing # index after which to insert the new node
        for i in 1:node
            for j in 1:n
                distance = distance_matrix[solution[i], j]
                if distance < smallest_distance
                    smallest_distance = distance
                    smallest_distance_node = j
                    from = i
                end
            end
        end
    end

    insert!()


end


distance_matrix = create_distance_matrix("lab1/TSPA.csv")
println("Running random search")
random_score, random_solution = find_best_solution(random_search, 200, distance_matrix)
println("best score: ", random_score)
println("solution: ", random_solution)

println("\nRunning nn_last search")
random_score, random_solution = find_best_solution(nn_last_node, 200, distance_matrix)
println("best score: ", random_score)
println("solution: ", random_solution)

println("\nRunning nn_all search")

