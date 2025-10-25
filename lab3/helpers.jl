using CSV, DataFrames, StatsBase, Plots, Random

function euclidean(a, b)
    return round(sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2))
end

function prepare_data(data_path)
    data = CSV.read(data_path, DataFrame; header=["x", "y", "w"])
    distance_matrix = zeros(nrow(data), nrow(data))
    costs = data[:, "w"]

    for i in 1:nrow(data)
        for j in 1:nrow(data)
            if i != j
                distance_matrix[i, j] = euclidean((data[i, "x"], data[i, "y"]), (data[j, "x"], data[j, "y"]))
            end
        end
    end

    return distance_matrix, costs
end

function calculate_cycle_length(solution, distance_matrix, costs)
    score = 0
    n = length(solution)
    for i in 1:n
        next_i = mod1(i + 1, n)  # Wraps around to 1 when i == n
        score += distance_matrix[solution[i], solution[next_i]] + costs[solution[next_i]]
    end
    return score
end

function random_start(distance_matrix, costs, weights, start)
    n = size(distance_matrix, 1)
    k = ceil(Int, n / 2)
    return randperm(n)[1:k]
end