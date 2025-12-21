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

function random_start(distance_matrix; rng=Random.default_rng())
    n = size(distance_matrix, 1)
    k = ceil(Int, n / 2)
    perm = randperm(rng, n)
    return perm[1:k]
end

function plot_best_solution(scores, solutions, data_path, title, save_path)
    best_solution_index = argmin(scores)
    best_solution = solutions[best_solution_index]

    data = CSV.read(data_path, DataFrame; header=["x", "y", "w"])
    all_costs = collect(data[:, "w"])
    min_cost, max_cost = extrema(all_costs)

    x_coords = [data[i, "x"] for i in best_solution]
    y_coords = [data[i, "y"] for i in best_solution]
    push!(x_coords, data[best_solution[1], "x"])
    push!(y_coords, data[best_solution[1], "y"])

    p = plot(x_coords, y_coords;
        line=:solid,
        linewidth=2,
        # marker=:circle,
        # markersize=6,
        label="Best Solution Path",
        xlabel="X",
        ylabel="Y",
        title=title,
        titlefont=font(8),
        legend=:best)

    scatter!(p,
        data[:, "x"],
        data[:, "y"];
        marker=:circle,
        markersize=4,
        marker_z=all_costs,
        color=cgrad([:blue, :red]),
        clims=(min_cost, max_cost),
        label="Node Cost",
        alpha=0.85,
        colorbar_title="Cost")

    mkpath(dirname(save_path))
    savefig(p, save_path)
end