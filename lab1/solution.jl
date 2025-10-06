using CSV, DataFrames, StatsBase

function euclidean(x1, y1, x2, y2)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

function create_distance_matrix(data_path)
    data = CSV.read("lab1/TSPA.csv", DataFrame; header=["x", "y", "w"])
    distance_matrix = zeros(nrow(data), nrow(data))

    for i in 1:nrow(data)
        for j in 1:nrow(data)
            if i != j
                x1 = data[i, "x"]
                y1 = data[i, "y"]
                x2 = data[j, "x"]
                y2 = data[j, "y"]
                distance_matrix[i, j] = euclidean(x1, y1, x2, y2) + data[i, "w"]
            end
        end
    end

    return distance_matrix
end

function find_best_solution(search_function, iterations, distance_matrix)
    best_score = Inf
    best_solution = nothing

    for i in 1:iterations
        score, solution = search_function(distance_matrix)
        if score < best_score
            best_score = score
            best_solution = solution
        end
    end

    return best_score, best_solution
end

function random_search(distance_matrix)
    n = ceil(Int, nrow(data) / 2)
    sampled_nodes = sample(1:size(distance_matrix)[1], n; replace=false)
    score = 0
    for i in 1:length(sampled_nodes)-1
        score += distance_matrix[sampled_nodes[i], sampled_nodes[i+1]]
    end

    score += distance_matrix[sampled_nodes[n], sampled_nodes[1]]
    return score, sampled_nodes
end


distance_matrix = create_distance_matrix("lab1/TSPA.csv")
println("Running random search")
random_score, random_solution = find_best_solution(random_search, 200, distance_matrix)
println("best score: ", random_score)
println("solution: ", random_solution)
