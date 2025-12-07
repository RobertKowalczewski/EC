using CSV
using DataFrames
using ProgressMeter

include("local-search.jl")
include("moves.jl")
include("helpers.jl")

println("TSPA")
data_path = "lab1/TSPA.csv"
distance_matrix, costs = prepare_data(data_path)

# solutions = []

# @showprogress 0.1 "runs" for i in 1:1000
#     x = random_start(distance_matrix)
#     obj, x = greedy_local_search(x, distance_matrix, costs, intra_two_edges_exchange, inter_two_nodes_exchange)
#     push!(x, obj)
#     push!(solutions, x)
# end

# println(solutions)

# # Convert the array of arrays to a DataFrame, where each solution is a row
# # The header will be automatically generated as x1, x2, ...
# df = DataFrame(permutedims(hcat(solutions...)), :auto)

# # Save the DataFrame to a CSV file
# CSV.write("solutions_TSPA.csv", df)

println("TSPB")
data_path = "lab1/TSPB.csv"
distance_matrix, costs = prepare_data(data_path)

# solutions = []

# @showprogress 0.1 "runs" for i in 1:1000
#     x = random_start(distance_matrix)
#     obj, x = greedy_local_search(x, distance_matrix, costs, intra_two_edges_exchange, inter_two_nodes_exchange)
#     push!(x, obj)
#     push!(solutions, x)
# end

# # Convert the array of arrays to a DataFrame, where each solution is a row
# # The header will be automatically generated as x1, x2, ...
# df = DataFrame(permutedims(hcat(solutions...)), :auto)

# # Save the DataFrame to a CSV file
# CSV.write("solutions_TSPB.csv", df)